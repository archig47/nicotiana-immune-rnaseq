# Complete Autonomous RNA-seq Agent Setup - Summary & Next Steps

## What You Now Have

You have a **complete, production-ready system** for generating and validating an RNA-seq Snakemake pipeline.

### Files Created (3 main documents):

1. **`autonomous-rna-seq-agent-prompt.md`** (500+ lines)
   - Complete system prompt for the Claude agent
   - Three-phase execution protocol (GENERATE → VALIDATE → OPTIMIZE)
   - Biological validation checklists
   - Parameter testing framework
   - Error handling protocols
   - **This is the specification the agent will follow**

2. **`master-invocation-guide.md`** (300+ lines)
   - Step-by-step instructions for invoking the agent
   - How to copy/paste the specification
   - How to extract generated files
   - Troubleshooting guide
   - **This is your user manual**

3. **`invoke-rna-seq-agent.sh`** (bash script)
   - Automated script to prepare the invocation
   - Copies invocation to clipboard
   - Platform-specific (macOS/Linux)
   - **This is a convenience tool**

---

## Your Next Steps (In Order)

### STEP 1: Prepare Project Directory (2 minutes)

```bash
# Navigate to your project
cd ~/rna-seq-pipeline-project

# Verify files are there
ls -la Snakefile config.yaml README.md 2>/dev/null || echo "Not yet created (that's OK)"

# Copy the agent documents here
cp /mnt/user-data/outputs/autonomous-rna-seq-agent-prompt.md ./
cp /mnt/user-data/outputs/master-invocation-guide.md ./
cp /mnt/user-data/outputs/invoke-rna-seq-agent.sh ./

# Make the invocation script executable
chmod +x invoke-rna-seq-agent.sh

# Verify
ls -la *.md invoke-rna-seq-agent.sh

# Commit
git add autonomous-rna-seq-agent-prompt.md master-invocation-guide.md invoke-rna-seq-agent.sh
git commit -m "Add autonomous agent system prompts and invocation guides"
```

### STEP 2: Review the Agent Specification (10 minutes)

```bash
# Read the complete agent specification
cat autonomous-rna-seq-agent-prompt.md

# Key sections to understand:
# - PHASE 1: GENERATION (what files will be created)
# - PHASE 2: VALIDATION (what checks will be performed)
# - PHASE 3: OPTIMIZATION (how parameters will be tested)
# - YOUR VALIDATION CHECKLIST (what "success" looks like)
```

### STEP 3: Invoke the Agent (2 minutes setup + 10-15 minutes execution)

**Option A: RECOMMENDED - Use the automated script**

```bash
cd ~/rna-seq-pipeline-project

# Run the invocation helper
bash invoke-rna-seq-agent.sh

# Script will:
# 1. Prepare the complete command
# 2. Tell you how to copy it
# 3. Instruct you to paste into Claude terminal
# 4. Save to /tmp/invocation_command.txt for easy access
```

**Option B: Manual invocation**

```bash
cd ~/rna-seq-pipeline-project

# Open Claude in your terminal
claude

# In Claude, paste the text below (all of it):
# ============================================

You are the Autonomous RNA-seq Snakemake Pipeline Generator & Validator.

PHASES:
1. GENERATE - Create complete Snakemake pipeline
2. VALIDATE - Test on Pombo 2019 data (SRP118889)
3. OPTIMIZE - Test parameter variations

COMPLETE SPECIFICATION:

[PASTE THE ENTIRE CONTENT OF: autonomous-rna-seq-agent-prompt.md]

Available tools and versions:
- fasterq-dump (SRA Toolkit) v3.2.1
- Falco v1.2.5
- FastQC v0.12.7
- MultiQC v1.33
- Trim Galore v0.6.11
- STAR v2.7.11b
- featureCounts (Subread) v2.1.1
- DESeq2 v1.50.2

BEGIN PHASE 1: GENERATION

Output each as delimited code block:
1. === SNAKEFILE START/END ===
2. === CONFIG.YAML START/END ===
3. === ENVIRONMENT.YML START/END ===
4. === README.MD START/END ===
5. === TEST_DATA.SH START/END ===
6. === VALIDATION_CHECKLIST.MD START/END ===
7. === GITIGNORE START/END ===

When complete, state: "PHASE 1 COMPLETE"

Then describe PHASE 2 and PHASE 3.

# ============================================
```

### STEP 4: Wait for Claude (10-15 minutes)

Claude will now:
- **Read** the complete specification (3-5 min)
- **Generate** the Snakefile (2-3 min)
- **Generate** config, environment, README files (2-3 min)
- **Output** validation checklists and descriptions (2-3 min)

**Total: ~10-15 minutes**

You'll see output like:
```
=== SNAKEFILE START ===
# Snakemake workflow for RNA-seq analysis
# [lots of lines]
=== SNAKEFILE END ===

=== CONFIG.YAML START ===
# Configuration for RNA-seq pipeline
# [lots of lines]
=== CONFIG.YAML END ===

[etc.]
```

### STEP 5: Extract Generated Files (10 minutes)

Once Claude finishes, you'll have ~7 code blocks to extract.

**For each code block:**

```bash
# Create the file
cat > ~/rna-seq-pipeline-project/FILENAME << 'EOF'
[PASTE THE CONTENT HERE - everything between START and END]
EOF

# Verify it was created
ls -lh ~/rna-seq-pipeline-project/FILENAME
```

**Files to extract:**
1. Snakefile (~200-300 lines)
2. config.yaml (~80-120 lines)
3. environment.yml (~25-35 lines)
4. README.md (~150-200 lines)
5. test_data.sh (~50-80 lines) - make executable: `chmod +x test_data.sh`
6. VALIDATION_CHECKLIST.md (~150-200 lines)
7. .gitignore (~20-30 lines)

### STEP 6: Verify All Files (2 minutes)

```bash
cd ~/rna-seq-pipeline-project

# Check all files exist
ls -la Snakefile config.yaml environment.yml README.md test_data.sh VALIDATION_CHECKLIST.md .gitignore

# Verify file sizes (rough check for completeness)
wc -l Snakefile config.yaml environment.yml README.md

# Expected:
# Snakefile: 200-400 lines
# config.yaml: 80-150 lines
# environment.yml: 25-40 lines
# README.md: 100-200 lines
```

### STEP 7: Commit to Git (2 minutes)

```bash
cd ~/rna-seq-pipeline-project

# Add all generated files
git add Snakefile config.yaml environment.yml README.md test_data.sh VALIDATION_CHECKLIST.md .gitignore

# Commit with message
git commit -m "Add agent-generated Snakemake pipeline (PHASE 1 COMPLETE)"

# Verify
git log --oneline -5
```

### STEP 8: Read Documentation (10 minutes)

```bash
cd ~/rna-seq-pipeline-project

# Read the README
cat README.md

# Read the validation checklist
cat VALIDATION_CHECKLIST.md

# Review the agent's Phase 2 and Phase 3 descriptions
# (Claude should have provided these after generating files)
cat AGENT_OUTPUT_PHASE2_3.md  # if you saved it
```

### STEP 9: Test Pipeline Setup (5 minutes)

```bash
cd ~/rna-seq-pipeline-project

# Create conda environment
mamba env create -f environment.yml
# OR
conda env create -f environment.yml

# This will take 5-10 minutes (installing all tools)

# Activate environment
conda activate rna-seq-pipeline

# Verify tools installed
snakemake --version
falco --version
star --version
```

### STEP 10: Ready for Testing (Whenever you want)

```bash
# When ready to test on Pombo data:
cd ~/rna-seq-pipeline-project

# Download data
./test_data.sh
# OR manually:
fasterq-dump --split-files SRR6676954 SRR6676955 SRR6676956 SRR6676957 SRR6676958 SRR6676959

# Run pipeline
snakemake results/deseq2/volcano_plot.pdf --cores 20 --use-conda

# This will run the complete pipeline on Pombo data
# Check progress against VALIDATION_CHECKLIST.md at each step
```

---

## Timeline Summary

| Task | Time | Status |
|------|------|--------|
| Prepare project directory | 2 min | Do now ✅ |
| Review agent specification | 10 min | Do now ✅ |
| Invoke agent (setup) | 2 min | Do now ✅ |
| Claude generation | 10-15 min | Happens automatically |
| Extract generated files | 10 min | Do after Claude finishes |
| Commit to git | 2 min | Do after extraction |
| Read documentation | 10 min | Do when you have time |
| Test pipeline setup | 5-10 min | Do before testing data |
| **TOTAL SETUP** | **~50-60 min** | |
| **Test on Pombo data** | **30-60 min** | (Next phase) |

---

## What Happens After STEP 10

Once you've tested the pipeline on Pombo data:

1. **Validate**: Go through VALIDATION_CHECKLIST.md item by item
2. **Troubleshoot**: Use ERROR HANDLING section if anything fails
3. **Optimize**: Run parameter variation tests (if needed)
4. **Scale up**: Test on your viral co-infection dataset (24-36 samples)
5. **Deploy**: Push to GitHub, document results, publish

---

## Key Files to Keep Track Of

```
~/rna-seq-pipeline-project/
├── Snakefile                        ← Generated by agent
├── config.yaml                      ← Generated by agent
├── environment.yml                  ← Generated by agent
├── README.md                        ← Generated by agent
├── test_data.sh                     ← Generated by agent
├── VALIDATION_CHECKLIST.md          ← Generated by agent
├── .gitignore                       ← Generated by agent
│
├── autonomous-rna-seq-agent-prompt.md    ← Agent spec (this session)
├── master-invocation-guide.md             ← User manual (this session)
├── invoke-rna-seq-agent.sh                ← Helper script (this session)
│
├── PIPELINE_SPECIFICATION.md        ← Earlier (saved copy)
├── PIPELINE_CONTEXT.txt             ← Earlier (your notes)
├── DECISIONS_LOG.md                 ← Earlier (your decisions)
│
├── data/
│   └── reads/                       ← Downloaded FASTQ files go here
├── results/                         ← Pipeline outputs go here
└── .git/                           ← Git repository
```

---

## What The Agent Will Check Biologically

The agent is programmed to validate:

✅ **Quality**: 
- All samples >10M reads
- GC content 45-48%
- Duplicates 37-64% (normal for RNA-seq)

✅ **Alignment**:
- Mapping rate >80%
- Strandedness confirmed (reverse-stranded)
- No batch effects in BAM files

✅ **Counting**:
- Sensible count distribution
- >95% assigned reads
- No obvious artifacts

✅ **DESeq2 Results**:
- PCA shows treatment separation
- DEGs include known immune genes (MAPK, WRKY, NLR, PR)
- log₂FC distribution symmetric
- No statistical artifacts

✅ **Biological Sense**:
- For Pst: Immune genes UP, housekeeping stable
- For viral data: Immune suppression patterns match expectations
- Replicates cluster together

---

## Questions Before You Start?

**Q: Can I run Claude in the background?**
A: Yes! Use `tmux` or `screen` to detach after pasting the invocation.

**Q: What if Claude times out?**
A: Ask it: "Continue generating the remaining files" in a new message.

**Q: Can I modify parameters after generation?**
A: YES! Edit `config.yaml` - all parameters are there, no need to regenerate.

**Q: What if validation fails on Pombo data?**
A: Check VALIDATION_CHECKLIST.md and ERROR HANDLING section for troubleshooting.

**Q: When do I test on my viral dataset?**
A: After successful validation on Pombo (6 samples). Then scale up to viral data (24-36 samples).

---

## You're Ready! 🚀

You now have a complete, autonomous system to:
1. **Generate** a production-grade Snakemake pipeline
2. **Validate** it on real data (Pombo 2019)
3. **Optimize** parameters based on results
4. **Deploy** to GitHub with full documentation

**Next action**: Run `invoke-rna-seq-agent.sh` or follow `master-invocation-guide.md`

Good luck! 🧬📊
