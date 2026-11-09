import os
import time
import json
import hashlib

import redis
from fastapi import FastAPI, Request, HTTPException
from pydantic import BaseModel, Field
from Bio.Seq import Seq

app = FastAPI(
    title="SeqLab API",
    description="DNA & Protein Sequence Analysis with Caching & Rate Limiting",
)

# ------------------------------------------------------------------ Redis setup
# Deployed (Render / Upstash): set REDIS_URL.
# Local / docker-compose: REDIS_HOST (defaults to localhost) is used instead.
REDIS_URL = os.getenv("REDIS_URL")
REDIS_HOST = os.getenv("REDIS_HOST", "localhost")

try:
    if REDIS_URL:
        redis_client = redis.from_url(
            REDIS_URL,
            decode_responses=True,
            socket_connect_timeout=2,
            socket_timeout=2,
        )
    else:
        redis_client = redis.Redis(
            host=REDIS_HOST,
            port=6379,
            db=0,
            decode_responses=True,
            socket_connect_timeout=2,
            socket_timeout=2,
        )
    redis_client.ping()
    REDIS_AVAILABLE = True
except Exception:
    redis_client = None
    REDIS_AVAILABLE = False

# ----------------------------------------------------------------------- Config
MAX_TOKENS = 5            # bucket capacity = max burst per client
REFILL_RATE = 1.0         # tokens added per second (sustained rate)
CACHE_TTL_SECONDS = 60
VALID_CHARS = set("ATGCUN")


class SequenceRequest(BaseModel):
    sequence: str = Field(..., description="DNA sequence string")


# ---------------------------------------------------------------------- Helpers
def validate_sequence(seq_str: str):
    if not seq_str:
        raise HTTPException(status_code=400, detail="Sequence cannot be empty.")
    if not all(c in VALID_CHARS for c in seq_str):
        raise HTTPException(
            status_code=400,
            detail="Invalid character found in sequence. Allowed: A, T, G, C, U, N.",
        )


def get_client_ip(request: Request) -> str:
    # Behind a proxy (Render) the real client IP is in X-Forwarded-For.
    forwarded = request.headers.get("x-forwarded-for")
    if forwarded:
        return forwarded.split(",")[0].strip()
    return request.client.host if request.client else "unknown"


def check_rate_limit(client_id: str) -> bool:
    """Token bucket stored in a Redis hash. Fails open if Redis is unavailable."""
    if not REDIS_AVAILABLE:
        return True
    try:
        key = f"rate_limit:{client_id}"
        now = time.time()
        bucket = redis_client.hgetall(key)
        if not bucket:
            tokens = MAX_TOKENS - 1
        else:
            elapsed = now - float(bucket["last_refill"])
            tokens = min(MAX_TOKENS, float(bucket["tokens"]) + elapsed * REFILL_RATE)
            if tokens < 1:
                return False
            tokens -= 1
        redis_client.hset(key, mapping={"tokens": tokens, "last_refill": now})
        redis_client.expire(key, 3600)
        return True
    except redis.RedisError:
        return True


def enforce_rate_limit(request: Request):
    if not check_rate_limit(get_client_ip(request)):
        raise HTTPException(
            status_code=429,
            detail="Rate limit exceeded.",
            headers={"Retry-After": "1"},  # one token refills per second
        )


def get_cache_key(operation: str, sequence: str) -> str:
    return hashlib.sha256(f"{operation}:{sequence.strip().upper()}".encode()).hexdigest()


def cache_get(key: str):
    """Single Redis GET (no exists()+get() race). Returns None on miss or Redis error."""
    if not REDIS_AVAILABLE:
        return None
    try:
        raw = redis_client.get(key)
        return json.loads(raw) if raw else None
    except (redis.RedisError, ValueError):
        return None


def cache_set(key: str, value: dict):
    if not REDIS_AVAILABLE:
        return
    try:
        redis_client.setex(key, CACHE_TTL_SECONDS, json.dumps(value))
    except redis.RedisError:
        pass


def run_cached(operation: str, request: Request, body: SequenceRequest, compute):
    """Shared flow: rate limit -> validate -> cache lookup -> compute -> cache store."""
    enforce_rate_limit(request)
    start = time.perf_counter()

    seq_str = body.sequence.strip().upper()
    validate_sequence(seq_str)

    key = get_cache_key(operation, seq_str)
    hit = cache_get(key)
    if hit is not None:
        return {
            "cached": True,
            "response_time_ms": round((time.perf_counter() - start) * 1000, 2),
            **hit,
        }

    try:
        result = compute(seq_str)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

    cache_set(key, result)
    return {
        "cached": False,
        "response_time_ms": round((time.perf_counter() - start) * 1000, 2),
        **result,
    }


# ------------------------------------------------------------ Analysis functions
def trim_to_codons(seq):
    return seq[: len(seq) - len(seq) % 3]


def compute_gc(seq_str: str) -> dict:
    seq = Seq(seq_str)
    gc = (seq.count("G") + seq.count("C")) / len(seq) * 100
    return {"sequence": seq_str, "length": len(seq), "gc_content_percent": round(gc, 2)}


def compute_reverse_complement(seq_str: str) -> dict:
    return {"original": seq_str, "reverse_complement": str(Seq(seq_str).reverse_complement())}


def compute_translation(seq_str: str) -> dict:
    return {"sequence": seq_str, "protein": str(trim_to_codons(Seq(seq_str)).translate())}


def compute_orfs(seq_str: str) -> dict:
    seq = Seq(seq_str)
    orfs = []
    for frame in range(3):
        translated = str(trim_to_codons(seq[frame:]).translate())
        for p in translated.split("*"):
            start = p.find("M")
            if start != -1:
                orfs.append(
                    {"frame": frame + 1, "protein_length": len(p) - start, "protein": p[start:]}
                )
    return {"sequence": seq_str, "orfs_found": len(orfs), "orfs": orfs}


# ----------------------------------------------------------------------- Routes
RESPONSES = {
    400: {"description": "Empty or invalid sequence"},
    429: {"description": "Rate limit exceeded"},
}


@app.get("/")
def home():
    return {"message": "Welcome to SeqLab!", "redis_connected": REDIS_AVAILABLE}


@app.post("/gc-content", responses=RESPONSES)
def calculate_gc(request: Request, body: SequenceRequest):
    return run_cached("gc-content", request, body, compute_gc)


@app.post("/reverse-complement", responses=RESPONSES)
def reverse_complement(request: Request, body: SequenceRequest):
    return run_cached("reverse-complement", request, body, compute_reverse_complement)


@app.post("/translate", responses=RESPONSES)
def translate_sequence(request: Request, body: SequenceRequest):
    return run_cached("translate", request, body, compute_translation)


@app.post("/orfs", responses=RESPONSES)
def find_orfs(request: Request, body: SequenceRequest):
    return run_cached("orfs", request, body, compute_orfs)