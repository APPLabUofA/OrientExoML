#!/usr/bin/env python3
"""
download_osf_eeg.py

Downloads the raw BrainVision EEG files (*_Orient_Exo_Exp.vhdr, .vmrk, .eeg)
for the included participants from OSF (https://osf.io/u6fgm/) using parallel
streams with size verification.
"""

import json
import os
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
EEG_DIR = PROJECT_ROOT / "data" / "eeg"
EEG_DIR.mkdir(parents=True, exist_ok=True)

INCLUDED_SUBS = {
    "004", "005", "007", "009", "010",
    "012", "013", "014", "015", "017",
    "018", "019", "021", "022", "024",
    "025", "026", "027", "029", "031",
    "032", "034", "035", "036", "037",
    "038", "039", "040"
}

OSF_EEG_FOLDER_URL = "https://api.osf.io/v2/nodes/u6fgm/files/osfstorage/62323377ac980a057d94fd77/?page[size]=100"


def fetch_osf_manifest():
    items = []
    url = OSF_EEG_FOLDER_URL
    while url:
        req = urllib.request.Request(url, headers={"User-Agent": "Python-OSF-Downloader"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            d = json.loads(resp.read().decode())
            items.extend(d.get("data", []))
            url = d.get("links", {}).get("next")
    return items


def download_file(name, size, url):
    dest = EEG_DIR / name
    if dest.exists() and dest.stat().st_size == size:
        return name, "cached", size

    tmp_dest = EEG_DIR / f"{name}.part"
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Python-OSF-Downloader"})
            with urllib.request.urlopen(req, timeout=120) as resp, open(tmp_dest, "wb") as out_f:
                while True:
                    chunk = resp.read(1024 * 1024)
                    if not chunk:
                        break
                    out_f.write(chunk)
            if tmp_dest.stat().st_size == size:
                tmp_dest.rename(dest)
                return name, "downloaded", size
        except Exception as e:
            time.sleep(2 * (attempt + 1))
            if attempt == 3:
                return name, f"failed: {e}", 0
    return name, "failed_size_mismatch", 0


def main():
    print("[download_osf_eeg] Fetching OSF file manifest...")
    items = fetch_osf_manifest()
    to_download = []
    for it in items:
        attr = it.get("attributes", {})
        name = attr.get("name", "")
        size = attr.get("size") or 0
        sub = name[:3]
        if sub in INCLUDED_SUBS and "_Exp." in name:
            dl_url = it["links"]["download"]
            to_download.append((name, size, dl_url))

    to_download.sort(key=lambda x: x[0])
    total_mb = sum(x[1] for x in to_download) / (1024 ** 2)
    print(f"[download_osf_eeg] Found {len(to_download)} files ({total_mb:.1f} MB) for included participants.")

    t0 = time.time()
    completed = 0
    bytes_done = 0
    with ThreadPoolExecutor(max_workers=6) as pool:
        futures = {pool.submit(download_file, n, s, u): (n, s) for n, s, u in to_download}
        for fut in as_completed(futures):
            name, status, sz = fut.result()
            completed += 1
            bytes_done += sz
            elapsed = time.time() - t0
            rate = (bytes_done / (1024 ** 2)) / max(elapsed, 0.1)
            print(f"  [{completed}/{len(to_download)}] {name}: {status} ({sz / (1024**2):.1f} MB) | {rate:.1f} MB/s", flush=True)

    print(f"[download_osf_eeg] Finished in {time.time() - t0:.1f}s.")


if __name__ == "__main__":
    main()

