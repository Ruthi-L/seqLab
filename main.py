import time
import hashlib
import json
import os
import redis
from fastapi import FastAPI, HTTPException, Request, status
from Bio.Seq import Seq



app = FastAPI(title="SeqLab API", description="DNA & Protein Sequence Analysis with Redis Caching")

# Connect to Redis (defaults to localhost:6379)
# decode_responses=True ensures Redis returns strings instead of bytes

# Get Redis host from environment variable, default to localhost for local testing
REDIS_HOST = os.getenv("REDIS_HOST", "localhost")

try:
    redis_client = redis.Redis(host=REDIS_HOST, port=6379, db=0, decode_responses=True)
    redis_client.ping()
    REDIS_AVAILABLE = True
except Exception:
    REDIS_AVAILABLE = False

def get_cache_key(operation: str, sequence: str) -> str:
    """Generate a unique hash key for a given sequence and operation."""
    combined = f"{operation}:{sequence.strip().upper()}"
    return hashlib.sha256(combined.encode()).hexdigest()

@app.get("/")
def home():
    return {
        "message": "Welcome to SeqLab! Go to /docs to test the API.",
        "redis_connected": REDIS_AVAILABLE
    }

# Rate limiter settings
MAX_TOKENS = 5       # Max requests allowed in burst
REFILL_RATE = 1.0    # 1 token added per second

def check_rate_limit(client_id: str) -> bool:
    """
    Hand-written Token-Bucket rate limiter using Redis.
    Returns True if request is allowed, False if rate-limited.
    """
    if not REDIS_AVAILABLE:
        return True # If Redis is down, don't block the user

    redis_key = f"rate_limit:{client_id}"
    now = time.time()

    # Get current bucket state from Redis (tokens and last update time)
    bucket = redis_client.hgetall(redis_key)
    
    if not bucket:
        # Brand new bucket: start full
        tokens = MAX_TOKENS - 1 # Take 1 for the current request
        last_refill = now
    else:
        tokens = float(bucket["tokens"])
        last_refill = float(bucket["last_refill"])

        # Calculate how many tokens to add based on elapsed time
        elapsed = now - last_refill
        tokens_to_add = elapsed * REFILL_RATE
        
        # Refill tokens, but don't exceed max capacity
        tokens = min(MAX_TOKENS, tokens + tokens_to_add)

        if tokens < 1:
            # Not enough tokens left! Block request.
            return False
        
        # Consume 1 token for this request
        tokens -= 1
        last_refill = now

    # Save updated bucket state back to Redis (expire after 1 hour of inactivity)
    redis_client.hset(redis_key, mapping={
        "tokens": tokens,
        "last_refill": last_refill
    })
    redis_client.expire(redis_key, 3600)
    
    return True


@app.post("/gc-content")
def calculate_gc(request: Request, sequence: str):
    """Calculate GC content with Rate Limiting and Redis Caching."""
    
    # 1. Check Rate Limit (using client IP address as identifier)
    client_ip = request.client.host if request.client else "unknown"
    if not check_rate_limit(client_ip):
        raise HTTPException(
            status_code=status.HTTP_429_TOO_MANY_REQUESTS,
            detail="Rate limit exceeded. Too many requests, please slow down."
        )

    start_time = time.time()
    seq_str = sequence.strip().upper()
    
    if not seq_str:
        raise HTTPException(status_code=400, detail="Sequence is empty.")

    cache_key = get_cache_key("gc-content", seq_str)
    
    # 2. Check Cache
    if REDIS_AVAILABLE and redis_client.exists(cache_key):
        cached_data = json.loads(redis_client.get(cache_key))
        response_time = round((time.time() - start_time) * 1000, 2)
        return {
            "cached": True,
            "response_time_ms": response_time,
            **cached_data
        }

    # 3. Compute if cache miss
    try:
        seq = Seq(seq_str)
        gc_count = seq.count("G") + seq.count("C")
        length = len(seq)
        gc_percentage = round((gc_count / length) * 100, 2)
        
        result = {
            "sequence": seq_str,
            "length": length,
            "gc_content_percent": gc_percentage
        }
        
        if REDIS_AVAILABLE:
            redis_client.setex(cache_key, 60, json.dumps(result))
            
        response_time = round((time.time() - start_time) * 1000, 2)
        return {
            "cached": False,
            "response_time_ms": response_time,
            **result
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))