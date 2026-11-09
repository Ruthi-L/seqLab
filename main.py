import os
import time
import hashlib
import json
import redis
from fastapi import FastAPI, Request, HTTPException, status
from pydantic import BaseModel, Field
from Bio.Seq import Seq

app = FastAPI(title="SeqLab API", description="DNA & Protein Sequence Analysis with Caching & Rate Limiting")

REDIS_HOST = os.getenv("REDIS_HOST", "localhost")

try:
    redis_client = redis.Redis(host=REDIS_HOST, port=6379, db=0, decode_responses=True)
    redis_client.ping()
    REDIS_AVAILABLE = True
except Exception:
    REDIS_AVAILABLE = False

MAX_TOKENS = 5
REFILL_RATE = 1.0

class SequenceRequest(BaseModel):
    sequence: str = Field(..., description="DNA sequence string")

def validate_sequence(seq_str: str):
    if not seq_str:
        raise HTTPException(status_code=400, detail="Sequence cannot be empty.")
    valid_chars = set("ATGCUN")
    if not all(c in valid_chars for c in seq_str):
        raise HTTPException(status_code=400, detail=f"Invalid character found in sequence. Allowed: A, T, G, C, U, N.")

def check_rate_limit(client_id: str) -> bool:
    if not REDIS_AVAILABLE:
        return True
    redis_key = f"rate_limit:{client_id}"
    now = time.time()
    bucket = redis_client.hgetall(redis_key)
    if not bucket:
        tokens = MAX_TOKENS - 1
        last_refill = now
    else:
        tokens = float(bucket["tokens"])
        last_refill = float(bucket["last_refill"])
        elapsed = now - last_refill
        tokens = min(MAX_TOKENS, tokens + (elapsed * REFILL_RATE))
        if tokens < 1:
            return False
        tokens -= 1
        last_refill = now
    redis_client.hset(redis_key, mapping={"tokens": tokens, "last_refill": last_refill})
    redis_client.expire(redis_key, 3600)
    return True

def get_cache_key(operation: str, sequence: str) -> str:
    return hashlib.sha256(f"{operation}:{sequence.strip().upper()}".encode()).hexdigest()

def enforce_rate_limit(request: Request):
    client_ip = request.client.host if request.client else "unknown"
    if not check_rate_limit(client_ip):
          raise HTTPException(status_code=429, detail="Rate limit exceeded.",
                      headers={"Retry-After": "5"})   # use your real refill time
      

def validate_sequence(seq_str: str):
    if not seq_str:
        raise HTTPException(status_code=400, detail="Sequence cannot be empty.")
    # Check for invalid non-nucleotide characters
    valid_chars = set("ATGCUN")
    if not all(c in valid_chars for c in seq_str):
        raise HTTPException(status_code=400, detail=f"Invalid character found in sequence. Allowed: A, T, G, C, U, N.")

@app.get("/")
def home():
    return {"message": "Welcome to SeqLab!", "redis_connected": REDIS_AVAILABLE}

@app.post("/gc-content")
def calculate_gc(request: Request, body: SequenceRequest):
    enforce_rate_limit(request)
    start_time = time.time()
    seq_str = body.sequence.strip().upper()
    validate_sequence(seq_str)
    
    cache_key = get_cache_key("gc-content", seq_str)
    if REDIS_AVAILABLE and redis_client.exists(cache_key):
        return {"cached": True, "response_time_ms": round((time.time() - start_time) * 1000, 2), **json.loads(redis_client.get(cache_key))}
    
    try:
        seq = Seq(seq_str)
        result = {"sequence": seq_str, "length": len(seq), "gc_content_percent": round(((seq.count("G") + seq.count("C")) / len(seq)) * 100, 2)}
        if REDIS_AVAILABLE:
            redis_client.setex(cache_key, 60, json.dumps(result))
        return {"cached": False, "response_time_ms": round((time.time() - start_time) * 1000, 2), **result}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/reverse-complement")
def reverse_complement(request: Request, body: SequenceRequest):
    enforce_rate_limit(request)
    start_time = time.time()
    seq_str = body.sequence.strip().upper()
    validate_sequence(seq_str)
    
    cache_key = get_cache_key("reverse-complement", seq_str)
    if REDIS_AVAILABLE and redis_client.exists(cache_key):
        return {"cached": True, "response_time_ms": round((time.time() - start_time) * 1000, 2), **json.loads(redis_client.get(cache_key))}
    
    try:
        seq = Seq(seq_str)
        result = {"original": seq_str, "reverse_complement": str(seq.reverse_complement())}
        if REDIS_AVAILABLE:
            redis_client.setex(cache_key, 60, json.dumps(result))
        return {"cached": False, "response_time_ms": round((time.time() - start_time) * 1000, 2), **result}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/translate")
def translate_sequence(request: Request, body: SequenceRequest):
    enforce_rate_limit(request)
    start_time = time.time()
    seq_str = body.sequence.strip().upper()
    validate_sequence(seq_str)
    
    cache_key = get_cache_key("translate", seq_str)
    if REDIS_AVAILABLE and redis_client.exists(cache_key):
        return {"cached": True, "response_time_ms": round((time.time() - start_time) * 1000, 2), **json.loads(redis_client.get(cache_key))}
    
    try:
        seq = Seq(seq_str)
        result = {"sequence": seq_str, "protein": str(seq.translate())}
        if REDIS_AVAILABLE:
            redis_client.setex(cache_key, 60, json.dumps(result))
        return {"cached": False, "response_time_ms": round((time.time() - start_time) * 1000, 2), **result}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/orfs")
def find_orfs(request: Request, body: SequenceRequest):
    enforce_rate_limit(request)
    start_time = time.time()
    seq_str = body.sequence.strip().upper()
    validate_sequence(seq_str)
    
    cache_key = get_cache_key("orfs", seq_str)
    if REDIS_AVAILABLE and redis_client.exists(cache_key):
        return {"cached": True, "response_time_ms": round((time.time() - start_time) * 1000, 2), **json.loads(redis_client.get(cache_key))}
    
    try:
        seq = Seq(seq_str)
        orfs = []
        for frame in range(3):
            translated = seq[frame:].translate(to_stop=False)
            for p in str(translated).split("*"):
                if "M" in p:
                    orfs.append({"frame": frame + 1, "protein_length": len(p[p.find("M"):]), "protein": p[p.find("M"):]})
        result = {"sequence": seq_str, "orfs_found": len(orfs), "orfs": orfs}
        if REDIS_AVAILABLE:
            redis_client.setex(cache_key, 60, json.dumps(result))
        return {"cached": False, "response_time_ms": round((time.time() - start_time) * 1000, 2), **result}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))