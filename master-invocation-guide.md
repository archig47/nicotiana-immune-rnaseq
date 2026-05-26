# Master Invocation Guide: Autonomous RNA-seq Agent

## TL;DR (If you just want to get started)

```bash
cd ~/rna-seq-pipeline-project

# Copy the entire content of autonomous-rna-seq-agent-prompt.md
cat autonomous-rna-seq-agent-prompt.md

# Open Claude in your terminal
claude

# Paste this (EXACTLY as shown):
```

**EXACT TEXT TO PASTE INTO CLAUDE:**

```
You are the Autonomous RNA-seq Snakemake Pipeline Generator & Validator.

PHASE 1: Generate complete, production-grade Snakemake pipeline
PHASE 2: Test on real data (Pombo 2019, SRP118889, 6 samples), validate at each step
PHASE 3: Optimize parameters, choose best based on biological results

COMPLETE SPECIFICATION BELOW:

=== SPECIFICATION START ===

[PASTE THE ENTIRE CONTENT OF: ~/rna-seq-pipeline-project/autonomous-rna-seq-agent-prompt.md]

=== SPECIFICATION END ===

TOOLS AVAILABLE:
- fasterq-dump (SRA Toolkit) v3.2.1
- Falco v1.2.5
- FastQC v0.12.7
- MultiQC v1.33
- Trim Galore v0.6.11
- STAR v2.7.11b
- featureCounts (Subread) v2.1.1
- DESeq2 v1.50.2

BEGIN PHASE 1: GENERATION

Generate and output ALL of the following as code blocks:

1. Snakefile (complete, with all rules)
2. config.yaml (with all parameters explained)
3. environment.yml (pinned versions)
4. README.md (complete usage guide)
5. test_data.sh (download script)
6. VALIDATION_CHECKLIST.md
7. .gitignore

After outputting files, state: "PHASE 1 COMPLETE"

Then describe PHASE 2 and PHASE 3 (don't actually run them, just describe what would happen)
```

---

## Step-by-Step Instructions

### Step 1: Prepare Your Project Directory

```bash
# If you haven't done this yet
mkdir -p ~/rna-seq-pipeline-project
cd ~/rna-seq-pipeline-project

# You should already have these files:
ls -la
# Should show:
# - PIPELINE_SPECIFICATION.md (from earlier)
# - autonomous-rna-seq-agent-prompt.md (just created)
# - PIPELINE_CONTEXT.txt (your context doc)
# - DECISIONS_LOG.md
# - .git/ (initialized)
```

### Step 2: Open Claude in Terminal

```bash
# Make sure you have Claude CLI installed
# If not: pip install anthropic

# Navigate to your project
cd ~/rna-seq-pipeline-project

# Open Claude
claude

# You should see a prompt like:
# "Type your message..."
```

### Step 3: Copy the Agent Specification

```bash
# In a separate terminal window (or before opening claude):

# Find the autonomous-rna-seq-agent-prompt.md file
cd ~/rna-seq-pipeline-project
cat autonomous-rna-seq-agent-prompt.md | wc -l

# It should be ~500+ lines (the complete specification)
# Copy the entire file to your clipboard
cat autonomous-rna-seq-agent-prompt.md > /tmp/agent_spec.txt

# (Now you have it saved in /tmp/agent_spec.txt for easy reference)
```

### Step 4: Paste the Invocation into Claude

In your Claude terminal session, paste:

```
You are the Autonomous RNA-seq Snakemake Pipeline Generator & Validator.

PHASES:
1. GENERATE: Complete Snakemake pipeline (all rules, configs, docs)
2. VALIDATE: Describe testing protocol on Pombo data
3. OPTIMIZE: Describe parameter optimization framework

You MUST FOLLOW the complete specification provided below EXACTLY.

=== SPECIFICATION START ===
```

Then paste the **entire content** of `autonomous-rna-seq-agent-prompt.md`

```
=== SPECIFICATION END ===

Available tools:
- fasterq-dump (SRA Toolkit) v3.2.1
- Falco v1.2.5
- FastQC v0.12.7
- MultiQC v1.33
- Trim Galore v0.6.11
- STAR v2.7.11b
- featureCounts (Subread) v2.1.1
- DESeq2 v1.50.2

BEGIN PHASE 1: GENERATION

Output all of the following as separate code blocks:

1. === SNAKEFILE START ===
   [Complete Snakefile content]
   === SNAKEFILE END ===

2. === CONFIG.YAML START ===
   [Complete config.yaml content]
   === CONFIG.YAML END ===

3. === ENVIRONMENT.YML START ===
   [Complete environment.yml content]
   === ENVIRONMENT.YML END ===

4. === README.MD START ===
   [Complete README.md content]
   === README.MD END ===

5. === TEST_DATA.SH START ===
   [Complete test_data.sh content]
   === TEST_DATA.SH END ===

6. === VALIDATION_CHECKLIST.MD START ===
   [Complete validation checklist]
   === VALIDATION_CHECKLIST.MD END ===

7. === GITIGNORE START ===
   [Complete .gitignore content]
   === GITIGNORE END ===

When complete, state exactly:
"PHASE 1 COMPLETE"

Then provide PHASE 2 and PHASE 3 descriptions.
```

---

### Step 5: Wait for Claude to Generate

Claude will now:
- Read the entire specification
- Generate all files
- Output them as code blocks
- Describe validation protocol
- Describe parameter optimization

This will take 5-15 minutes depending on Claude's response time.

---

### Step 6: Extract the Generated Files

Once Claude finishes, you'll see code blocks like:

```
=== SNAKEFILE START ===
# Snakemake workflow for RNA-seq analysis
# ... [content] ...
=== SNAKEFILE END ===
```

**For each file block:**

1. **Copy the content** (everything between START and END)
2. **Create the local file**:

```bash
# Example for Snakefile
cat > ~/rna-seq-pipeline-project/Snakefile << 'EOF'
[PASTE THE CONTENT HERE]
EOF

# Example for config.yaml
cat > ~/rna-seq-pipeline-project/config.yaml << 'EOF'
[PASTE THE CONTENT HERE]
EOF

# Repeat for all files:
# - config.yaml
# - environment.yml
# - README.md
# - test_data.sh (make executable: chmod +x test_data.sh)
# - VALIDATION_CHECKLIST.md
# - .gitignore
```

---

### Step 7: Verify Files Were Created

```bash
cd ~/rna-seq-pipeline-project

# Check that all files exist
ls -la Snakefile config.yaml environment.yml README.md test_data.sh VALIDATION_CHECKLIST.md .gitignore

# Should show:
# -rw-r--r--  Snakefile
# -rw-r--r--  config.yaml
# -rw-r--r--  environment.yml
# -rw-r--r--  README.md
# -rwxr-xr-x  test_data.sh  (executable)
# -rw-r--r--  VALIDATION_CHECKLIST.md
# -rw-r--r--  .gitignore

# Make test_data.sh executable if not already
chmod +x test_data.sh

# Verify file sizes (rough check)
wc -l Snakefile config.yaml environment.yml README.md
# Snakefile should be ~200-400 lines
# config.yaml should be ~50-100 lines
# environment.yml should be ~20-30 lines
# README.md should be ~100-150 lines
```

---

### Step 8: Commit to Git

```bash
cd ~/rna-seq-pipeline-project

git add Snakefile config.yaml environment.yml README.md test_data.sh VALIDATION_CHECKLIST.md .gitignore

git commit -m "Add agent-generated RNA-seq Snakemake pipeline"

git log --oneline
# Should show your new commit at the top
```

---

### Step 9: Read Claude's Validation & Optimization Descriptions

Claude should have also provided:

- **PHASE 2 VALIDATION PROTOCOL**: Step-by-step instructions for testing on Pombo data
- **PHASE 3 OPTIMIZATION FRAMEWORK**: How to test parameter variants

**Save these descriptions:**

```bash
# Create a file to store them
cat > ~/rna-seq-pipeline-project/AGENT_OUTPUT_PHASE2_3.md << 'EOF'
[PASTE CLAUDE'S PHASE 2 and PHASE 3 DESCRIPTIONS HERE]
EOF
```

---

### Step 10: Next Steps (When You're Ready)

```bash
cd ~/rna-seq-pipeline-project

# Read the README
cat README.md

# Read the validation checklist
cat VALIDATION_CHECKLIST.md

# When you're ready to test on Pombo data:
# (Follow the PHASE 2 instructions from Claude)

# Quick test: Dry run
snakemake -n --cores 1

# Or follow test_data.sh instructions:
./test_data.sh
```

---

## If Something Goes Wrong

### Claude stopped mid-generation?
Ask it: "Continue from where you left off" or "Generate the remaining files (config.yaml, environment.yml, etc.)"

### Files are incomplete?
Ask Claude: "Check the Snakefile for completeness. Are all rules defined? Does it have validation checkpoints?"

### Validation checklist is empty?
Ask Claude: "Generate a detailed validation checklist with PASS/FAIL criteria for each step"

### Need clarification on parameters?
Ask Claude: "Explain the Trim Galore parameters and why each was chosen" (etc. for any tool)

---

## What You Should Have at the End

✅ **Snakefile** - Ready to run (`snakemake --cores 20`)
✅ **config.yaml** - Ready to customize for your datasets
✅ **environment.yml** - Reproducible conda environment
✅ **README.md** - Complete documentation
✅ **test_data.sh** - Download and test on Pombo data
✅ **VALIDATION_CHECKLIST.md** - Know what success looks like
✅ **.gitignore** - Clean GitHub repo
✅ **PHASE 2 & 3 descriptions** - Know how to validate and optimize

---

## Quick Checklist

- [ ] Created ~/rna-seq-pipeline-project/
- [ ] Have autonomous-rna-seq-agent-prompt.md in project directory
- [ ] Opened claude in terminal
- [ ] Pasted the agent invocation (including full specification)
- [ ] Waited for Claude to complete
- [ ] Copied all generated files into local files
- [ ] Made test_data.sh executable (`chmod +x test_data.sh`)
- [ ] Committed to git (`git add -A && git commit`)
- [ ] Read README.md and VALIDATION_CHECKLIST.md
- [ ] Ready to test on Pombo data

---

## Timeline

| Action | Time |
|--------|------|
| Setup project directory | 2 min |
| Open Claude | 1 min |
| Paste specification | 2 min |
| Claude generates files | 10-15 min |
| Extract files | 5-10 min |
| Commit to git | 1 min |
| **TOTAL** | **~20-30 min** |

Then you'll be ready to test the pipeline!

---

## Questions?

**Q: Can I run this in the background?**
A: Yes - open Claude in a `tmux` or `screen` session, paste the prompt, then detach. Claude will keep generating.

**Q: How long does Phase 1 take?**
A: Usually 10-15 minutes for Claude to think through and generate all files.

**Q: What if I need to modify parameters later?**
A: Edit `config.yaml` - that's what it's for! No need to regenerate.

**Q: Can I test this on my viral dataset right away?**
A: No - test on Pombo first (6 samples, public data). Then expand to viral dataset (24-36 samples).

**Q: What if validation fails?**
A: Check VALIDATION_CHECKLIST.md and ERROR HANDLING section in agent prompt for troubleshooting.

---

You're ready! 🚀
