"""
Provenance tracking for ML pipeline reproducibility
Generates sidecar .meta.json files for each output
"""
import json
import hashlib
import os
import datetime
import sys
import subprocess
from pathlib import Path


def file_sha256(path: str) -> str:
    """Calculate SHA256 hash of a file"""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def get_git_commit() -> str:
    """Get current git commit hash"""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=10
        )
        return result.stdout.strip()
    except Exception:
        return "unknown"


def get_git_branch() -> str:
    """Get current git branch"""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            capture_output=True,
            text=True,
            timeout=10
        )
        return result.stdout.strip()
    except Exception:
        return "unknown"


def write_provenance(
    output_path: str,
    input_paths: list,
    config: dict,
    rule_name: str,
    script_path: str
) -> None:
    """Write provenance metadata sidecar file"""
    
    meta = {
        "output_file": os.path.basename(output_path),
        "output_path": os.path.abspath(output_path),
        "output_hash": file_sha256(output_path),
        "rule_name": rule_name,
        "script_path": script_path,
        "created_utc": datetime.datetime.utcnow().isoformat() + "Z",
        "runtime": {
            "python_version": sys.version,
            "platform": sys.platform
        },
        "git": {
            "commit": get_git_commit(),
            "branch": get_git_branch()
        },
        "config": {
            "disease": config.get("project", {}).get("disease", "Cancer"),
            "organism": config.get("project", {}).get("organism", "9606"),
            "string_version": config.get("string", {}).get("version_tag", "v12.0"),
            "rwr_restart_prob": config.get("propagation", {}).get("rwr", {}).get("restart_prob", 0.7),
            "rwr_top_k": config.get("propagation", {}).get("rwr", {}).get("top_k", 500),
            "clustering_resolution": config.get("clustering", {}).get("resolution", 1.0)
        },
        "inputs": []
    }
    
    # Track input files
    for inp_path in input_paths:
        if os.path.exists(inp_path):
            meta["inputs"].append({
                "path": inp_path,
                "hash": file_sha256(inp_path)
            })
    
    # Write sidecar file
    meta_path = output_path + ".meta.json"
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, sort_keys=True)
    
    return meta_path
