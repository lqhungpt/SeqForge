#!/usr/bin/env python3
"""
Viral Assembly HTML Report Generator
Usage: python3 generate_report.py --prefix PREFIX --outdir OUTDIR --depth DEPTH_TXT [options]
"""

import argparse
import os
import json
import subprocess
from datetime import datetime

def parse_flagstat(path):
    stats = {}
    if not os.path.exists(path):
        return stats
    with open(path) as f:
        for line in f:
            line = line.strip()
            if "in total" in line or ("total" in line and "QC" in line):
                stats["total"] = line.split()[0]
            elif "mapped" in line and "%" in line and "mate" not in line:
                parts = line.split()
                stats["mapped"] = parts[0]
                stats["mapped_pct"] = parts[4].strip("(")
            elif "properly paired" in line:
                parts = line.split()
                stats["properly_paired"] = parts[0]
                stats["properly_paired_pct"] = parts[4].strip("(")
    return stats

def parse_fastp_json(path):
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)

def read_depth(path, downsample=3000):
    positions, depths = [], []
    if not os.path.exists(path):
        return positions, depths
    with open(path) as f:
        lines = f.readlines()
    step = max(1, len(lines) // downsample)
    for i in range(0, len(lines), step):
        parts = lines[i].strip().split()
        if len(parts) >= 3:
            positions.append(int(parts[1]))
            depths.append(int(parts[2]))
    return positions, depths

def read_fasta(path):
    """Read FASTA and return header + sequence string."""
    if not os.path.exists(path):
        return "", ""
    header, seq = "", []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line
            else:
                seq.append(line)
    return header, "".join(seq)

def parse_contigs(path):
    """Return list of (name, length) for each contig in FASTA."""
    contigs = []
    if not os.path.exists(path):
        return contigs
    name, seq = "", []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    contigs.append((name, len("".join(seq))))
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name:
        contigs.append((name, len("".join(seq))))
    return contigs

def parse_fastqc_data(fastqc_dir, prefix_hint=""):
    """Parse FastQC per-base quality and per-sequence quality from fastqc_data.txt files."""
    results = []
    if not os.path.exists(fastqc_dir):
        return results
    for root, dirs, files in os.walk(fastqc_dir):
        for fname in files:
            if fname == "fastqc_data.txt":
                path = os.path.join(root, fname)
                sample = os.path.basename(root).replace("_fastqc", "")
                pb_qual, ps_qual = [], []
                section = None
                with open(path) as f:
                    for line in f:
                        line = line.rstrip()
                        if line.startswith(">>Per base sequence quality"):
                            section = "pb"
                        elif line.startswith(">>Per sequence quality scores"):
                            section = "ps"
                        elif line.startswith(">>END_MODULE"):
                            section = None
                        elif section == "pb" and not line.startswith("#"):
                            parts = line.split("\t")
                            if len(parts) >= 2:
                                pb_qual.append({"pos": parts[0], "mean": float(parts[1])})
                        elif section == "ps" and not line.startswith("#"):
                            parts = line.split("\t")
                            if len(parts) >= 2:
                                try:
                                    ps_qual.append({"q": int(float(parts[0])), "count": float(parts[1])})
                                except ValueError:
                                    pass
                results.append({"sample": sample, "pb_qual": pb_qual, "ps_qual": ps_qual})
    return results

def get_tool_version(cmd):
    """Run a version command and return first line of output."""
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=10)
        out = (result.stdout + result.stderr).strip()
        return out.split("\n")[0] if out else "N/A"
    except Exception:
        return "N/A"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix",           required=True)
    parser.add_argument("--outdir",           required=True)
    parser.add_argument("--depth",            required=True)
    parser.add_argument("--consensus_len",    default="0")
    parser.add_argument("--consensus_gc",     default="0")
    parser.add_argument("--consensus_n",      default="0")
    parser.add_argument("--avg_depth",        default="0")
    parser.add_argument("--min_depth",        default="0")
    parser.add_argument("--max_depth",        default="0")
    parser.add_argument("--covered_pct",      default="0")
    parser.add_argument("--low_depth_bases",  default="0")
    parser.add_argument("--n_count",          default="0")
    parser.add_argument("--contig_count",     default="0")
    parser.add_argument("--snps_corrected",   default="0")
    parser.add_argument("--gaps_filled",      default="0")
    parser.add_argument("--final_mapped",     default="0")
    parser.add_argument("--elapsed_min",      default="0")
    parser.add_argument("--elapsed_sec",      default="0")
    parser.add_argument("--threads",          default="0")
    parser.add_argument("--denovo_reads",      default="0")
    parser.add_argument("--fastqc_raw_dir",    default="")
    parser.add_argument("--fastqc_final_dir",  default="")
    args = parser.parse_args()

    p   = args.prefix
    out = args.outdir

    # Paths
    fastp_json    = os.path.join(out, "01_trimmed",           f"{p}_fastp.json")
    ref_flagstat  = os.path.join(out, "02_mapped",            f"{p}_ref_flagstat.txt")
    scaf_flagstat = os.path.join(out, "05_ragtag_scaffold",   f"{p}_scaffold_flagstat.txt")
    fin_flagstat  = os.path.join(out, "07_coverage_map",      f"{p}_final_flagstat.txt")
    consensus_fa  = os.path.join(out, "06_final_consensus",   f"{p}_final_consensus.fasta")
    contigs_fa    = os.path.join(out, "04_shovill_denovo",    f"{p}_contigs.fasta")
    fastqc_raw_dir   = args.fastqc_raw_dir   or os.path.join(out, "00_fastqc_raw")
    fastqc_final_dir = args.fastqc_final_dir or os.path.join(out, "08_fastqc_final")

    # Data
    positions, depths = read_depth(args.depth)
    fastp_data   = parse_fastp_json(fastp_json)
    ref_flag     = parse_flagstat(ref_flagstat)
    scaf_flag    = parse_flagstat(scaf_flagstat)
    fin_flag     = parse_flagstat(fin_flagstat)
    fa_header, fa_seq = read_fasta(consensus_fa)
    contigs_list     = parse_contigs(contigs_fa)
    fastqc_raw       = parse_fastqc_data(fastqc_raw_dir)
    fastqc_final     = parse_fastqc_data(fastqc_final_dir)

    # FASTP stats
    bf  = fastp_data.get("summary", {}).get("before_filtering", {})
    af  = fastp_data.get("summary", {}).get("after_filtering",  {})
    total_before = int(bf.get("total_reads", 0))
    total_after  = int(af.get("total_reads", 0))
    q30_before   = f"{float(bf.get('q30_rate', 0))*100:.1f}" if bf else "N/A"
    q30_after    = f"{float(af.get('q30_rate', 0))*100:.1f}" if af else "N/A"
    gc_before    = f"{float(bf.get('gc_content', 0))*100:.1f}" if bf else "N/A"
    gc_after     = f"{float(af.get('gc_content', 0))*100:.1f}" if af else "N/A"

    # Software versions
    versions = {
        "fastp":         get_tool_version("fastp --version 2>&1 | head -1"),
        "bwa-mem2":      get_tool_version("bwa-mem2 version 2>&1 | head -1"),
        "samtools":      get_tool_version("samtools version 2>&1 | head -1"),
        "SPAdes":        get_tool_version("spades.py --version 2>&1 | head -1"),
        "Shovill":       get_tool_version("shovill --version 2>&1 | head -1"),
        "RagTag":        get_tool_version("ragtag.py --version 2>&1 | head -1"),
        "Pilon":         get_tool_version("pilon --version 2>&1 | head -1"),
        "minimap2":      get_tool_version("minimap2 --version 2>&1 | head -1"),
        "FastQC":        get_tool_version("fastqc --version 2>&1 | head -1"),
        "Python":        get_tool_version("python3 --version 2>&1 | head -1"),
    }

    # Chart data
    depth_labels = json.dumps(positions)
    depth_values = json.dumps(depths)
    avg_line     = [float(args.avg_depth)] * len(positions)
    avg_values   = json.dumps(avg_line)
    gc_val       = float(args.consensus_gc)

    # Runtime
    elapsed_min = args.elapsed_min
    elapsed_sec = args.elapsed_sec
    threads     = args.threads

    # Software version table rows
    version_rows = "\n".join(
        f'<tr><td>{tool}</td><td><code>{ver}</code></td></tr>'
        for tool, ver in versions.items()
    )

    # FastQC chart data builder
    def fqc_pb_data(fqc_list):
        """Per-base quality: positions and mean Q scores."""
        if not fqc_list:
            return "[]", "[]"
        d = fqc_list[0]["pb_qual"]
        labels = json.dumps([x["pos"] for x in d])
        values = json.dumps([x["mean"] for x in d])
        return labels, values

    def fqc_ps_data(fqc_list):
        """Per-sequence quality: Q score bins and counts."""
        if not fqc_list:
            return "[]", "[]"
        d = fqc_list[0]["ps_qual"]
        labels = json.dumps([x["q"] for x in d])
        values = json.dumps([x["count"] for x in d])
        return labels, values

    raw_pb_labels,   raw_pb_values   = fqc_pb_data(fastqc_raw)
    raw_ps_labels,   raw_ps_values   = fqc_ps_data(fastqc_raw)
    final_pb_labels, final_pb_values = fqc_pb_data(fastqc_final)
    final_ps_labels, final_ps_values = fqc_ps_data(fastqc_final)

    # Q score breakdown (Q<20, Q20-29, Q30+)
    def q_breakdown(fqc_list):
        if not fqc_list:
            return 0, 0, 0
        d = fqc_list[0]["ps_qual"]
        low = sum(x["count"] for x in d if x["q"] < 20)
        mid = sum(x["count"] for x in d if 20 <= x["q"] < 30)
        high = sum(x["count"] for x in d if x["q"] >= 30)
        return int(low), int(mid), int(high)

    raw_q_low, raw_q_mid, raw_q_high       = q_breakdown(fastqc_raw)
    final_q_low, final_q_mid, final_q_high = q_breakdown(fastqc_final)

    # Contig size table rows
    contig_rows = ""
    total_contig_len = 0
    for i, (cname, clen) in enumerate(contigs_list, 1):
        total_contig_len += clen
        contig_rows += f"<tr><td>Contig {i}</td><td>{cname}</td><td>{clen:,} bp</td></tr>\n"
    if not contig_rows:
        contig_rows = "<tr><td colspan=\'3\'>No contig data available</td></tr>"

    # Contig bar chart data
    contig_names  = [f"Contig {i}" for i in range(1, len(contigs_list)+1)]
    contig_sizes  = [c[1] for c in contigs_list]

    # Sequence display (wrap at 60 chars)
    seq_wrapped = "\n".join(fa_seq[i:i+60] for i in range(0, len(fa_seq), 60)) if fa_seq else "N/A"
    seq_display = f"{fa_header}\n{seq_wrapped}"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Viral Assembly Report — {p}</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/4.4.1/chart.umd.min.js"></script>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         background: #f0f4f8; color: #2d3748; }}
  header {{ background: linear-gradient(135deg, #1a365d 0%, #2b6cb0 100%);
            color: white; padding: 28px 40px; }}
  header h1 {{ font-size: 24px; font-weight: 700; }}
  header p  {{ font-size: 13px; opacity: .8; margin-top: 4px; }}
  .container {{ max-width: 1200px; margin: 0 auto; padding: 28px 20px; }}
  .grid-2 {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin-bottom: 20px; }}
  .grid-3 {{ display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 20px; margin-bottom: 20px; }}
  .grid-4 {{ display: grid; grid-template-columns: repeat(4,1fr); gap: 16px; margin-bottom: 20px; }}
  .grid-5 {{ display: grid; grid-template-columns: repeat(5,1fr); gap: 16px; margin-bottom: 20px; }}
  .card {{ background: white; border-radius: 12px; padding: 20px;
           box-shadow: 0 1px 3px rgba(0,0,0,.08); }}
  .card h2 {{ font-size: 14px; font-weight: 600; color: #718096;
              text-transform: uppercase; letter-spacing: .05em; margin-bottom: 14px; }}
  .stat-card {{ background: white; border-radius: 12px; padding: 18px;
                box-shadow: 0 1px 3px rgba(0,0,0,.08); text-align: center; }}
  .stat-card .val {{ font-size: 28px; font-weight: 700; color: #2b6cb0; }}
  .stat-card .lbl {{ font-size: 12px; color: #718096; margin-top: 4px; }}
  .stat-card.green .val {{ color: #276749; }}
  .stat-card.red .val   {{ color: #c53030; }}
  .stat-card.orange .val{{ color: #c05621; }}
  .stat-card.purple .val{{ color: #553c9a; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 13px; }}
  th {{ background: #ebf4ff; padding: 9px 12px; text-align: left;
        font-weight: 600; color: #2b6cb0; border-bottom: 2px solid #bee3f8; }}
  td {{ padding: 8px 12px; border-bottom: 1px solid #e2e8f0; }}
  tr:last-child td {{ border-bottom: none; }}
  tr:hover td {{ background: #f7fafc; }}
  .badge {{ display: inline-block; padding: 2px 8px; border-radius: 9999px;
            font-size: 11px; font-weight: 600; }}
  .badge.good {{ background: #c6f6d5; color: #276749; }}
  .badge.warn {{ background: #fefcbf; color: #744210; }}
  .badge.bad  {{ background: #fed7d7; color: #c53030; }}
  .chart-wrap {{ position: relative; height: 300px; }}
  .chart-wrap-sm {{ position: relative; height: 200px; }}
  .section-title {{ font-size: 18px; font-weight: 700; color: #1a365d;
                    margin: 28px 0 12px; padding-left: 12px;
                    border-left: 4px solid #2b6cb0; }}
  .seq-box {{ font-family: "Courier New", monospace; font-size: 11px;
              background: #1a202c; color: #68d391; padding: 16px;
              border-radius: 8px; white-space: pre-wrap; word-break: break-all;
              max-height: 300px; overflow-y: auto; line-height: 1.6; }}
  .copy-btn {{ display: inline-block; margin-top: 10px; padding: 8px 18px;
               background: #2b6cb0; color: white; border: none;
               border-radius: 6px; cursor: pointer; font-size: 13px; font-weight: 600; }}
  .copy-btn:hover {{ background: #2c5282; }}
  .blast-btn {{ display: inline-block; margin-top: 10px; margin-left: 8px;
                padding: 8px 18px; background: #276749; color: white;
                border: none; border-radius: 6px; cursor: pointer;
                font-size: 13px; font-weight: 600; text-decoration: none; }}
  .blast-btn:hover {{ background: #22543d; }}
  .runtime-grid {{ display: grid; grid-template-columns: repeat(3,1fr); gap: 16px; margin-bottom: 20px; }}
  code {{ background: #edf2f7; padding: 1px 5px; border-radius: 3px;
          font-family: "Courier New", monospace; font-size: 12px; }}
  footer {{ text-align: center; padding: 24px; font-size: 12px; color: #a0aec0; }}
  @media print {{
    .copy-btn, .blast-btn {{ display: none; }}
    .seq-box {{ max-height: none; }}
  }}
</style>
</head>
<body>

<header>
  <h1>🧬 Viral Assembly Report — {p}</h1>
  <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} &nbsp;|&nbsp;
     Pipeline: FASTP → bwa-mem2 → Shovill/SPAdes → RagTag → Pilon → Coverage</p>
</header>

<div class="container">

  <!-- KEY METRICS -->
  <p class="section-title">Key Metrics</p>
  <div class="grid-5">
    <div class="stat-card green">
      <div class="val">{int(args.consensus_len):,}</div>
      <div class="lbl">Consensus Length (bp)</div>
    </div>
    <div class="stat-card {'green' if float(args.avg_depth) >= 100 else 'orange'}">
      <div class="val">{args.avg_depth}×</div>
      <div class="lbl">Average Coverage Depth</div>
    </div>
    <div class="stat-card {'green' if args.consensus_n == '0' else 'red'}">
      <div class="val">{args.consensus_n}</div>
      <div class="lbl">Remaining Ns</div>
    </div>
    <div class="stat-card purple">
      <div class="val">{args.contig_count}</div>
      <div class="lbl">De Novo Contigs</div>
    </div>
  </div>

  <!-- RUNTIME -->
  <p class="section-title">Run Performance</p>
  <div class="runtime-grid">
    <div class="stat-card">
      <div class="val">{elapsed_min}m {elapsed_sec}s</div>
      <div class="lbl">Total Runtime</div>
    </div>
    <div class="stat-card">
      <div class="val">{threads}</div>
      <div class="lbl">CPU Threads Used</div>
    </div>
    <div class="stat-card">
      <div class="val">{datetime.now().strftime('%Y-%m-%d')}</div>
      <div class="lbl">Run Date</div>
    </div>
  </div>

  <!-- CONSENSUS STATS -->
  <p class="section-title">Consensus Sequence</p>
  <div class="grid-2">
    <div class="card">
      <h2>GC Content</h2>
      <div class="chart-wrap-sm"><canvas id="gcChart"></canvas></div>
    </div>
    <div class="card">
      <h2>Assembly Statistics</h2>
      <table>
        <tr><th>Metric</th><th>Value</th></tr>
        <tr><td>Total length</td><td>{int(args.consensus_len):,} bp</td></tr>
        <tr><td>GC content</td><td>{args.consensus_gc} %</td></tr>
        <tr><td>Remaining Ns</td><td>{args.consensus_n}</td></tr>
        <tr><td>Contigs (de novo)</td><td>{args.contig_count}</td></tr>
        <tr><td>SNPs corrected (Pilon)</td><td>{args.snps_corrected}</td></tr>
        <tr><td>Gaps filled (Pilon)</td><td>{args.gaps_filled}</td></tr>
      </table>
    </div>
  </div>

  <!-- COVERAGE DEPTH -->
  <p class="section-title">Coverage Depth Profile</p>
  <div class="card" style="margin-bottom:20px">
    <h2>Per-position Sequencing Depth across Final Consensus</h2>
    <div class="chart-wrap"><canvas id="depthChart"></canvas></div>
  </div>

  <!-- COVERAGE STATS -->
  <p class="section-title">Coverage Statistics</p>
  <div class="grid-2" style="margin-bottom:20px">
    <div class="card">
      <h2>Depth Summary</h2>
      <table>
        <tr><th>Metric</th><th>Value</th><th>Status</th></tr>
        <tr><td>Average depth</td><td>{args.avg_depth} ×</td>
            <td><span class="badge {'good' if float(args.avg_depth)>=100 else 'warn'}">
            {'High' if float(args.avg_depth)>=100 else 'Low'}</span></td></tr>
        <tr><td>Minimum depth</td><td>{args.min_depth} ×</td>
            <td><span class="badge {'good' if float(args.min_depth)>=10 else 'warn'}">
            {'OK' if float(args.min_depth)>=10 else 'Low'}</span></td></tr>
        <tr><td>Maximum depth</td><td>{args.max_depth} ×</td>
            <td><span class="badge good">OK</span></td></tr>
        <tr><td>Coverage breadth</td><td>{args.covered_pct} %</td>
            <td><span class="badge {'good' if float(args.covered_pct)>=99 else 'warn'}">
            {'Complete' if float(args.covered_pct)>=99 else 'Partial'}</span></td></tr>
        <tr><td>Bases &lt; 10× depth</td><td>{args.low_depth_bases} bp</td>
            <td><span class="badge {'good' if int(args.low_depth_bases)==0 else 'warn'}">
            {'None' if int(args.low_depth_bases)==0 else 'Present'}</span></td></tr>
      </table>
    </div>
    <div class="card">
      <h2>Depth Distribution</h2>
      <div class="chart-wrap-sm"><canvas id="depthDistChart"></canvas></div>
    </div>
  </div>

  <!-- READ QC -->
  <p class="section-title">Read Quality (FASTP)</p>
  <div class="grid-2" style="margin-bottom:20px">
    <div class="card">
      <h2>Read Statistics</h2>
      <table>
        <tr><th>Metric</th><th>Before Trimming</th><th>After Trimming</th></tr>
        <tr><td>Total reads</td><td>{total_before:,}</td><td>{total_after:,}</td></tr>
        <tr><td>Q30 rate (%)</td><td>{q30_before}</td><td>{q30_after}</td></tr>
        <tr><td>GC content (%)</td><td>{gc_before}</td><td>{gc_after}</td></tr>
        <tr><td>Reads removed</td><td colspan="2">{total_before - total_after:,} ({((total_before-total_after)/total_before*100) if total_before else 0:.1f}%)</td></tr>
      </table>
    </div>
    <div class="card">
      <h2>Read Retention</h2>
      <div class="chart-wrap-sm"><canvas id="readChart"></canvas></div>
    </div>
  </div>

  <!-- MAPPING STATS -->
  <p class="section-title">Mapping Statistics</p>
  <div class="card" style="margin-bottom:20px">
    <h2>Reads Mapped at Each Pipeline Stage</h2>
    <table>
      <tr><th>Stage</th><th>Total Reads</th><th>Mapped</th><th>Mapped %</th><th>Properly Paired</th><th>Status</th></tr>
      <tr>
        <td>Step 2 · Reference mapping</td>
        <td>{ref_flag.get('total','N/A')}</td>
        <td>{ref_flag.get('mapped','N/A')}</td>
        <td>{ref_flag.get('mapped_pct','N/A')}</td>
        <td>{ref_flag.get('properly_paired','N/A')}</td>
        <td><span class="badge good">Complete</span></td>
      </tr>
      <tr>
        <td>Step 6 · Scaffold mapping (Pilon)</td>
        <td>{scaf_flag.get('total','N/A')}</td>
        <td>{scaf_flag.get('mapped','N/A')}</td>
        <td>{scaf_flag.get('mapped_pct','N/A')}</td>
        <td>{scaf_flag.get('properly_paired','N/A')}</td>
        <td><span class="badge good">Complete</span></td>
      </tr>
      <tr>
        <td>Step 7 · Final consensus mapping</td>
        <td>{fin_flag.get('total','N/A')}</td>
        <td>{fin_flag.get('mapped','N/A')}</td>
        <td>{fin_flag.get('mapped_pct','N/A')}</td>
        <td>{fin_flag.get('properly_paired','N/A')}</td>
        <td><span class="badge good">Complete</span></td>
      </tr>
    </table>
  </div>

  <!-- FASTQC RAW READS -->
  <p class="section-title">FastQC — Raw Reads Quality</p>
  <div class="grid-2" style="margin-bottom:20px">
    <div class="card">
      <h2>Per-base Sequence Quality (Raw)</h2>
      <div class="chart-wrap"><canvas id="rawPbChart"></canvas></div>
    </div>
    <div class="card">
      <h2>Per-sequence Quality Score Distribution (Raw)</h2>
      <div class="chart-wrap"><canvas id="rawPsChart"></canvas></div>
    </div>
  </div>
  <div class="card" style="margin-bottom:20px">
    <h2>Q Score Breakdown — Raw Reads</h2>
    <table>
      <tr><th>Q Score Range</th><th>Illumina Category</th><th>Read Count</th><th>Interpretation</th></tr>
      <tr><td>&lt; Q20</td><td><span class="badge bad">Low quality</span></td><td>{raw_q_low:,}</td><td>Error rate &gt; 1% — usually discarded by trimming</td></tr>
      <tr><td>Q20 – Q29</td><td><span class="badge warn">Acceptable</span></td><td>{raw_q_mid:,}</td><td>Error rate 0.1–1% — trimmed reads may fall here</td></tr>
      <tr><td>≥ Q30</td><td><span class="badge good">High quality</span></td><td>{raw_q_high:,}</td><td>Error rate &lt; 0.1% — ideal for assembly</td></tr>
    </table>
  </div>

  <!-- FASTQC FINAL MAPPED READS -->
  <p class="section-title">FastQC — Final Mapped Reads Quality</p>
  <div class="grid-2" style="margin-bottom:20px">
    <div class="card">
      <h2>Per-base Sequence Quality (Final BAM)</h2>
      <div class="chart-wrap"><canvas id="finalPbChart"></canvas></div>
    </div>
    <div class="card">
      <h2>Per-sequence Quality Score Distribution (Final BAM)</h2>
      <div class="chart-wrap"><canvas id="finalPsChart"></canvas></div>
    </div>
  </div>
  <div class="card" style="margin-bottom:20px">
    <h2>Q Score Breakdown — Final Mapped Reads</h2>
    <table>
      <tr><th>Q Score Range</th><th>Illumina Category</th><th>Read Count</th><th>Interpretation</th></tr>
      <tr><td>&lt; Q20</td><td><span class="badge bad">Low quality</span></td><td>{final_q_low:,}</td><td>Error rate &gt; 1%</td></tr>
      <tr><td>Q20 – Q29</td><td><span class="badge warn">Acceptable</span></td><td>{final_q_mid:,}</td><td>Error rate 0.1–1%</td></tr>
      <tr><td>≥ Q30</td><td><span class="badge good">High quality</span></td><td>{final_q_high:,}</td><td>Error rate &lt; 0.1%</td></tr>
    </table>
  </div>

  <!-- CONTIG SIZES -->
  <p class="section-title">De Novo Contig Sizes (Shovill/SPAdes)</p>
  <div class="grid-2" style="margin-bottom:20px">
    <div class="card">
      <h2>Contig Summary</h2>
      <table>
        <tr><th>#</th><th>Contig Name</th><th>Length</th></tr>
        {contig_rows}
      </table>
    </div>
    <div class="card">
      <h2>Contig Size Comparison</h2>
      <div class="chart-wrap-sm"><canvas id="contigChart"></canvas></div>
    </div>
  </div>

  <!-- FINAL CONSENSUS SEQUENCE -->
  <p class="section-title">Final Consensus Sequence</p>
  <div class="card" style="margin-bottom:20px">
    <h2>FASTA sequence — copy to BLAST at NCBI</h2>
    <button class="copy-btn" onclick="copySeq()">📋 Copy FASTA</button>
    <a class="blast-btn"
       href="https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?PAGE_TYPE=BlastSearch&PROGRAM=blastn&PAGE=Nucleotides"
       target="_blank">🔗 Open NCBI BLASTn</a>
    <div class="seq-box" id="seqBox">{seq_display}</div>
  </div>

  <!-- SOFTWARE VERSIONS -->
  <p class="section-title">Software Versions</p>
  <div class="card" style="margin-bottom:20px">
    <h2>All tools used in this pipeline</h2>
    <table>
      <tr><th>Tool</th><th>Version / Info</th></tr>
      {version_rows}
    </table>
  </div>

</div>

<footer>
  Viral Assembly Pipeline &nbsp;|&nbsp; Sample: {p} &nbsp;|&nbsp;
  Runtime: {elapsed_min}m {elapsed_sec}s &nbsp;|&nbsp;
  Threads: {threads} &nbsp;|&nbsp;
  {datetime.now().strftime('%Y-%m-%d')}
</footer>

<script>
// ── Copy sequence ─────────────────────────────────────────────────────────
function copySeq() {{
  const text = document.getElementById('seqBox').innerText;
  navigator.clipboard.writeText(text).then(() => {{
    const btn = document.querySelector('.copy-btn');
    btn.textContent = '✅ Copied!';
    setTimeout(() => btn.textContent = '📋 Copy FASTA', 2000);
  }});
}}

// ── Depth line chart ──────────────────────────────────────────────────────
new Chart(document.getElementById('depthChart'), {{
  type: 'line',
  data: {{
    labels: {depth_labels},
    datasets: [
      {{ label: 'Coverage depth', data: {depth_values},
         borderColor: '#2b6cb0', backgroundColor: 'rgba(43,108,176,.1)',
         borderWidth: 1, pointRadius: 0, fill: true, tension: 0.2 }},
      {{ label: 'Average ({args.avg_depth}×)', data: {avg_values},
         borderColor: '#e53e3e', borderWidth: 1.5,
         borderDash: [6,3], pointRadius: 0, fill: false }}
    ]
  }},
  options: {{
    responsive: true, maintainAspectRatio: false,
    plugins: {{ legend: {{ position: 'top' }} }},
    scales: {{
      x: {{ title: {{ display: true, text: 'Genomic Position (bp)' }} }},
      y: {{ title: {{ display: true, text: 'Depth (×)' }}, beginAtZero: true }}
    }}
  }}
}});

// ── GC doughnut ───────────────────────────────────────────────────────────
new Chart(document.getElementById('gcChart'), {{
  type: 'doughnut',
  data: {{
    labels: ['GC ({gc_val}%)', 'AT ({100-gc_val:.2f}%)'],
    datasets: [{{ data: [{gc_val}, {100-gc_val}],
                  backgroundColor: ['#2b6cb0','#bee3f8'], borderWidth: 0 }}]
  }},
  options: {{
    responsive: true, maintainAspectRatio: false,
    plugins: {{ legend: {{ position: 'bottom' }} }}, cutout: '70%'
  }}
}});


// ── Read retention + mapping bar ─────────────────────────────────────────
new Chart(document.getElementById('readChart'), {{
  type: 'bar',
  data: {{
    labels: ['Before trimming', 'After trimming', 'Mapped to ref', 'To de novo assembly', 'Mapped to consensus'],
    datasets: [{{ label: 'Reads',
                  data: [{total_before}, {total_after},
                         {ref_flag.get('mapped', 0) or 0},
                         {args.denovo_reads},
                         {fin_flag.get('mapped', 0) or 0}],
                  backgroundColor: ['#bee3f8','#2b6cb0','#9ae6b4','#f6ad55','#276749'],
                  borderRadius: 6 }}]
  }},
  options: {{
    responsive: true, maintainAspectRatio: false,
    plugins: {{ legend: {{ display: false }} }},
    scales: {{ y: {{ beginAtZero: true,
                     title: {{ display: true, text: 'Read count' }} }} }}
  }}
}});

// ── Depth distribution histogram ──────────────────────────────────────────
const rawDepths = {depth_values};
const maxD = Math.min({args.max_depth}, 5000);
const binSize = Math.max(1, Math.round(maxD / 20));
const bins = Array(20).fill(0);
rawDepths.forEach(d => {{
  const i = Math.min(Math.floor(d / binSize), 19);
  bins[i]++;
}});
const binLabels = bins.map((_, i) => (i * binSize) + '–' + ((i+1) * binSize));
new Chart(document.getElementById('depthDistChart'), {{
  type: 'bar',
  data: {{
    labels: binLabels,
    datasets: [{{ label: 'Windows', data: bins,
                  backgroundColor: '#2b6cb0', borderRadius: 3, borderWidth: 0 }}]
  }},
  options: {{
    responsive: true, maintainAspectRatio: false,
    plugins: {{ legend: {{ display: false }} }},
    scales: {{
      x: {{ title: {{ display: true, text: 'Depth (×)' }},
             ticks: {{ maxRotation: 45 }} }},
      y: {{ title: {{ display: true, text: 'Windows' }}, beginAtZero: true }}
    }}
  }}
}});
// ── FastQC per-base quality helper ───────────────────────────────────────
function qColor(q) {{
  if (q >= 28) return '#276749';
  if (q >= 20) return '#d69e2e';
  return '#c53030';
}}
function pbColors(vals) {{ return vals.map(v => qColor(v)); }}

// ── Raw per-base quality ──────────────────────────────────────────────────
const rawPbLabels = {raw_pb_labels};
const rawPbValues = {raw_pb_values};
if (rawPbLabels.length > 0) {{
  new Chart(document.getElementById('rawPbChart'), {{
    type: 'bar',
    data: {{
      labels: rawPbLabels,
      datasets: [{{ label: 'Mean Q score', data: rawPbValues,
                    backgroundColor: pbColors(rawPbValues), borderWidth: 0 }}]
    }},
    options: {{
      responsive: true, maintainAspectRatio: false,
      plugins: {{ legend: {{ display: false }},
                  annotation: {{ annotations: {{
                    q30: {{ type: 'line', yMin: 30, yMax: 30,
                            borderColor: '#2b6cb0', borderWidth: 1.5,
                            borderDash: [4,4], label: {{ content: 'Q30', display: true }} }},
                    q20: {{ type: 'line', yMin: 20, yMax: 20,
                            borderColor: '#e53e3e', borderWidth: 1.5,
                            borderDash: [4,4], label: {{ content: 'Q20', display: true }} }}
                  }} }} }},
      scales: {{
        x: {{ title: {{ display: true, text: 'Position in read (bp)' }},
               ticks: {{ maxTicksLimit: 20 }} }},
        y: {{ min: 0, max: 42, title: {{ display: true, text: 'Mean Quality (Phred)' }} }}
      }}
    }}
  }});
}} else {{
  document.getElementById('rawPbChart').parentElement.innerHTML =
    '<p style="color:var(--color-text-secondary);padding:20px">FastQC raw data not available — run pipeline with raw reads</p>';
}}

// ── Raw per-sequence quality ──────────────────────────────────────────────
const rawPsLabels = {raw_ps_labels};
const rawPsValues = {raw_ps_values};
if (rawPsLabels.length > 0) {{
  new Chart(document.getElementById('rawPsChart'), {{
    type: 'line',
    data: {{
      labels: rawPsLabels,
      datasets: [{{ label: 'Read count', data: rawPsValues,
                    borderColor: '#2b6cb0', backgroundColor: 'rgba(43,108,176,.1)',
                    borderWidth: 2, pointRadius: 2, fill: true }}]
    }},
    options: {{
      responsive: true, maintainAspectRatio: false,
      plugins: {{ legend: {{ display: false }} }},
      scales: {{
        x: {{ title: {{ display: true, text: 'Phred Quality Score' }} }},
        y: {{ beginAtZero: true, title: {{ display: true, text: 'Read Count' }} }}
      }}
    }}
  }});
}} else {{
  document.getElementById('rawPsChart').parentElement.innerHTML =
    '<p style="color:var(--color-text-secondary);padding:20px">FastQC raw data not available</p>';
}}

// ── Final per-base quality ────────────────────────────────────────────────
const finalPbLabels = {final_pb_labels};
const finalPbValues = {final_pb_values};
if (finalPbLabels.length > 0) {{
  new Chart(document.getElementById('finalPbChart'), {{
    type: 'bar',
    data: {{
      labels: finalPbLabels,
      datasets: [{{ label: 'Mean Q score', data: finalPbValues,
                    backgroundColor: pbColors(finalPbValues), borderWidth: 0 }}]
    }},
    options: {{
      responsive: true, maintainAspectRatio: false,
      plugins: {{ legend: {{ display: false }} }},
      scales: {{
        x: {{ title: {{ display: true, text: 'Position in read (bp)' }},
               ticks: {{ maxTicksLimit: 20 }} }},
        y: {{ min: 0, max: 42, title: {{ display: true, text: 'Mean Quality (Phred)' }} }}
      }}
    }}
  }});
}} else {{
  document.getElementById('finalPbChart').parentElement.innerHTML =
    '<p style="color:var(--color-text-secondary);padding:20px">FastQC final data not available</p>';
}}

// ── Final per-sequence quality ────────────────────────────────────────────
const finalPsLabels = {final_ps_labels};
const finalPsValues = {final_ps_values};
if (finalPsLabels.length > 0) {{
  new Chart(document.getElementById('finalPsChart'), {{
    type: 'line',
    data: {{
      labels: finalPsLabels,
      datasets: [{{ label: 'Read count', data: finalPsValues,
                    borderColor: '#276749', backgroundColor: 'rgba(39,103,73,.1)',
                    borderWidth: 2, pointRadius: 2, fill: true }}]
    }},
    options: {{
      responsive: true, maintainAspectRatio: false,
      plugins: {{ legend: {{ display: false }} }},
      scales: {{
        x: {{ title: {{ display: true, text: 'Phred Quality Score' }} }},
        y: {{ beginAtZero: true, title: {{ display: true, text: 'Read Count' }} }}
      }}
    }}
  }});
}} else {{
  document.getElementById('finalPsChart').parentElement.innerHTML =
    '<p style="color:var(--color-text-secondary);padding:20px">FastQC final data not available</p>';
}}

// ── Contig size bar chart ─────────────────────────────────────────────────
new Chart(document.getElementById('contigChart'), {{
  type: 'bar',
  data: {{
    labels: {json.dumps(contig_names)},
    datasets: [{{ label: 'Length (bp)',
                  data: {json.dumps(contig_sizes)},
                  backgroundColor: ['#2b6cb0','#68d391','#f6ad55','#fc8181'],
                  borderRadius: 6 }}]
  }},
  options: {{
    responsive: true, maintainAspectRatio: false,
    indexAxis: 'y',
    plugins: {{ legend: {{ display: false }} }},
    scales: {{
      x: {{ beginAtZero: true,
             title: {{ display: true, text: 'Length (bp)' }} }}
    }}
  }}
}});
</script>
</body>
</html>"""

    out_path = os.path.join(args.outdir, f"{args.prefix}_assembly_report.html")
    with open(out_path, "w") as f:
        f.write(html)
    print(f"Report written: {out_path}")

if __name__ == "__main__":
    main()
