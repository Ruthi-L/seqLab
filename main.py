from fastapi import FastAPI, HTTPException
from Bio.Seq import Seq

app = FastAPI(title="SeqLab API", description="DNA & Protein Sequence Analysis")

@app.get("/")
def home():
    return {"message": "Welcome to SeqLab! Go to /docs to test the API."}

@app.post("/gc-content")
def calculate_gc(sequence: str):
    """Calculate the GC content percentage of a DNA sequence."""
    try:
        seq = Seq(sequence.strip().upper())
        gc_count = seq.count("G") + seq.count("C")
        length = len(seq)
        if length == 0:
            raise ValueError("Sequence is empty.")
        gc_percentage = (gc_count / length) * 100
        return {"sequence": str(seq), "length": length, "gc_content_percent": round(gc_percentage, 2)}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/reverse-complement")
def reverse_complement(sequence: str):
    """Get the reverse complement of a DNA sequence."""
    try:
        seq = Seq(sequence.strip().upper())
        # Biopython makes reverse complementing effortless
        rev_comp = seq.reverse_complement()
        return {"original": str(seq), "reverse_complement": str(rev_comp)}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/translate")
def translate_sequence(sequence: str):
    """Translate a DNA sequence into an amino acid sequence (protein)."""
    try:
        seq = Seq(sequence.strip().upper())
        # Translate DNA to protein. stop_symbol="*" stops at the first stop codon by default if desired, 
        # but let's translate the whole thing standardly.
        protein = seq.translate()
        return {"sequence": str(seq), "protein": str(protein)}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/orf-finder")
def find_orfs(sequence: str):
    """Find Open Reading Frames (ORFs) in a DNA sequence."""
    try:
        seq = Seq(sequence.strip().upper())
        orfs = []
        
        # Search all 3 forward reading frames
        for strand, n in [(seq, +1)]:
            for frame in range(3):
                length = len(strand)
                # Translate starting from the frame offset
                translated = strand[frame:].translate(to_stop=False)
                # Find positions of potential protein chunks between start (M) and stop (*)
                protein_seqs = str(translated).split("*")
                current_pos = frame
                for p in protein_seqs:
                    if "M" in p:
                        # Find where 'M' starts in this chunk
                        m_idx = p.find("M")
                        orf_protein = p[m_idx:]
                        orfs.append({
                            "frame": frame + 1,
                            "protein_length": len(orf_protein),
                            "protein": orf_protein
                        })
                    current_pos += len(p) + 1 # +1 for the stop codon
                    
        return {"sequence": str(seq), "orfs_found": len(orfs), "orfs": orfs}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))