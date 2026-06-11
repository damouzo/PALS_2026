# WORKSHOP PLAN: Beyond the Data: Workflow Management and GRNs in scMultiome
Format: Masterclass + Live Execution Demo (No student laptops required)
Duration: 60 Minutes

## 1. Class Structure & Timing

### Block 1: The Biological Framework & The Endgame (15 Mins)
- Concept: Introduction to scMultiome networks (scRNA-seq + scATAC-seq).
- The Problem: How do we assign directional weights to a structural chromatin map? 
- The Tool: CellOracle. Explain the goal: moving from a static network diagram to a dynamic system capable of "In Silico Perturbations" (simulating a gene knockout and modeling cell-fate shifts via vector fields).
- Objective: Show them the *final biological plot* first, so they know exactly what problem we are solving computationally.

### Block 2: The Infrastructure Crisis & Nextflow Architecture (15 Mins)
- The Legacy Code: Open the terminal/editor and show them the realistic, messy setup: 
  - `GRN_dataProcess.R` (Seurat/Signac/Cicero)
  - `GRN_analysis.py` (Scanpy/CellOracle)
  - `GRN.sh` (The Bash wrapper attempting custom YAML parsing hacks).
- Point out the vulnerabilities: What happens if Python fails at minute 59? What if paths change? How do we handle conflicting R and Python package dependencies cleanly?
- Introduce Nextflow Concepts: Dataflow paradigm, processes, channels, polyglot tasks (R and Python in the same workflow), and the magic of `-resume`.

### Block 3: The Live-Demo Symphony (30 Mins)
- Action 1 (The Clean Code): Show the refactored DSL2 Nextflow pipeline. Point out how clean the pipeline logic looks when separated from the data scripts.
- Action 2 (First Execution): Launch `nextflow run main.nf -profile test` live on your laptop. Show the progress bars tracking the execution of the R preprocessing step followed by the Python network inference.
- Action 3 (The WOW Factor): Change a key Python parameter on the fly (e.g., network link pruning percentile from 98 to 95). Re-run the pipeline using the `-resume` flag. 
- The Lesson: Show them how Nextflow instantly caches and skips the heavy R/Cicero preprocessing step, running only the modified Python step in seconds.
- Action 4 (Output Inspection): Open the generated directory, show the clean logs, and display the final simulated cell-fate plots. Explain that this entire reproducible framework will be shared with them via GitHub after the talk.