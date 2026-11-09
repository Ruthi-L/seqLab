import random, time, statistics

from fastapi import requests

URL = "http://localhost:8000/orfs"
N = 20          # number of trials
LENGTH = 6000   # bases per random sequence
PAUSE = 1.0     # seconds between requests (raise this if you still hit 429s)

def post(seq):
    while True:
        r = requests.post(URL, json={"sequence": seq})
        if r.status_code == 429:      # rate limited: wait for bucket refill, retry
            time.sleep(5)
            continue
        r.raise_for_status()
        return r.json()

def rand_seq():
    return "".join(random.choices("ACGT", k=LENGTH))

post(rand_seq())      # warm-up: first call has one-time overhead, so discard it

miss, hit = [], []
for _ in range(N):
    seq = rand_seq()                  # fresh sequence, so guaranteed cache miss
    a = post(seq); time.sleep(PAUSE)
    b = post(seq); time.sleep(PAUSE)  # identical sequence, so cache hit
    assert a["cached"] is False and b["cached"] is True, "cache behaviour unexpected"
    miss.append(a["response_time_ms"])
    hit.append(b["response_time_ms"])

print(f"Sequence length: {LENGTH} bases, runs: {N}")
print(f"Cache miss: median {statistics.median(miss):.2f} ms, mean {statistics.mean(miss):.2f} ms")
print(f"Cache hit : median {statistics.median(hit):.2f} ms, mean {statistics.mean(hit):.2f} ms")
print(f"Speedup (median): {statistics.median(miss) / statistics.median(hit):.1f}x")