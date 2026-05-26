#!/bin/bash

# invoke-rna-seq-agent.sh
# This script prepares and displays the complete invocation command for the RNA-seq agent
# Usage: bash invoke-rna-seq-agent.sh

set -e

PROJECT_DIR="${1:-.}"
SPEC_FILE="$PROJECT_DIR/autonomous-rna-seq-agent-prompt.md"

# Check if specification file exists
if [ ! -f "$SPEC_FILE" ]; then
    echo "❌ ERROR: Could not find $SPEC_FILE"
    echo "Make sure you're in the project directory with the agent specification."
    exit 1
fi

clear

cat << 'BANNER'
╔════════════════════════════════════════════════════════════════════════════╗
║                                                                            ║
║        RNA-seq Snakemake Pipeline - Autonomous Agent Invocation           ║
║                                                                            ║
║  This script prepares the complete invocation for the Claude agent.       ║
║  The agent will:                                                          ║
║    1. GENERATE: Complete Snakemake pipeline (Snakefile, configs, etc.)   ║
║    2. VALIDATE: Describe testing on Pombo 2019 data (SRP118889)          ║
║    3. OPTIMIZE: Describe parameter optimization framework                ║
║                                                                            ║
╚════════════════════════════════════════════════════════════════════════════╝

BANNER

echo ""
echo "📋 STEP 1: Opening Claude Terminal Session"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Please open Claude in your terminal by running:"
echo ""
echo "  $ claude"
echo ""
read -p "Press ENTER when Claude is ready and showing a prompt..."

echo ""
echo "✅ Great! Claude terminal is ready."
echo ""
echo "📋 STEP 2: Preparing the Invocation Command"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Count lines in spec
SPEC_LINES=$(wc -l < "$SPEC_FILE")
echo "Agent specification size: $SPEC_LINES lines"
echo ""

cat << 'INSTRUCTIONS'
You will now paste a LONG command into Claude. The command includes:
  - The complete agent system prompt
  - The full RNA-seq pipeline specification
  - Tool versions and availability
  - Invocation instructions

TOTAL PASTE LENGTH: ~800-1000 lines

This is NORMAL. Claude can handle this.

INSTRUCTIONS

echo ""
echo "📋 STEP 3: Copy the Invocation Command"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Create a temporary file with the full invocation
TEMP_INVOCATION=$(mktemp)

cat > "$TEMP_INVOCATION" << 'INVOCATION_START'
You are the Autonomous RNA-seq Snakemake Pipeline Generator & Validator.

YOUR ROLE:
PHASE 1: GENERATE - Create complete, production-grade Snakemake pipeline
PHASE 2: VALIDATE - Describe testing protocol on real data (Pombo 2019)
PHASE 3: OPTIMIZE - Describe parameter optimization framework

FOLLOW THE COMPLETE SPECIFICATION BELOW EXACTLY.

You MUST output all generated files as clearly delimited code blocks.

=== FULL SPECIFICATION START ===

INVOCATION_START

cat "$SPEC_FILE" >> "$TEMP_INVOCATION"

cat >> "$TEMP_INVOCATION" << 'INVOCATION_END'

=== FULL SPECIFICATION END ===

AVAILABLE TOOLS AND VERSIONS:
  - fasterq-dump (SRA Toolkit) v3.2.1
  - Falco v1.2.5
  - FastQC v0.12.7
  - MultiQC v1.33
  - Trim Galore v0.6.11
  - STAR v2.7.11b
  - featureCounts (Subread) v2.1.1
  - DESeq2 v1.50.2 (R v4.5.3)

NOW EXECUTE PHASE 1: GENERATION

Output each file as a clearly delimited code block:

1. === SNAKEFILE START ===
   [Complete Snakefile - all rules included]
   === SNAKEFILE END ===

2. === CONFIG.YAML START ===
   [Complete config.yaml with all parameters explained]
   === CONFIG.YAML END ===

3. === ENVIRONMENT.YML START ===
   [Complete environment.yml with pinned versions]
   === ENVIRONMENT.YML END ===

4. === README.MD START ===
   [Comprehensive README.md]
   === README.MD END ===

5. === TEST_DATA.SH START ===
   [Complete test_data.sh script]
   === TEST_DATA.SH END ===

6. === VALIDATION_CHECKLIST.MD START ===
   [Validation checklist from specification]
   === VALIDATION_CHECKLIST.MD END ===

7. === GITIGNORE START ===
   [Standard .gitignore for RNA-seq pipeline]
   === GITIGNORE END ===

REQUIREMENTS FOR EACH FILE:
  - Snakefile: Include ALL rules (QC, alignment, counting, DESeq2, plotting)
  - Snakefile: Include validation checkpoints at each major stage
  - config.yaml: Comment each parameter with justification
  - environment.yml: Pin all tool versions exactly as specified
  - README.md: Include quick start, troubleshooting, biological background
  - test_data.sh: Download Pombo data, verify pipeline runs
  - VALIDATION_CHECKLIST.md: PASS/FAIL criteria for each step
  - .gitignore: Results/, logs/, .snakemake/, conda environment

When PHASE 1 is complete, state exactly:
"PHASE 1 COMPLETE"

Then describe PHASE 2 (VALIDATION):
  - Step-by-step testing on Pombo data
  - Validation checklist for each stage
  - How to interpret results biologically
  - Error handling and troubleshooting

Then describe PHASE 3 (OPTIMIZATION):
  - Parameter testing framework
  - How to compare options
  - Biological criteria for choosing best parameters
  - Expected metrics at each optimization level

BEGIN NOW.

INVOCATION_END

# Display instructions
cat << 'COPY_INSTRUCTIONS'

📋 STEP 3: Copy & Paste into Claude
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

The complete invocation command is ready. It's quite long (~800 lines).

OPTION A (RECOMMENDED): Copy from file to clipboard
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

On macOS:
  cat /tmp/invocation_command.txt | pbcopy

On Linux (with xclip):
  cat /tmp/invocation_command.txt | xclip -selection clipboard

On Linux (with xsel):
  cat /tmp/invocation_command.txt | xsel --clipboard

Then:
  - Switch to Claude terminal
  - Paste: Ctrl+Shift+V (or Cmd+V on macOS)
  - Press ENTER


OPTION B: View in editor and manual copy
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  vim /tmp/invocation_command.txt
  (or: nano, less, cat, or your preferred editor)

Then:
  - Select all (Ctrl+A in vim/nano)
  - Copy (Ctrl+C)
  - Paste into Claude terminal


OPTION C: File location for reference
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Complete invocation saved to:
  /tmp/invocation_command.txt

You can reference it later with:
  cat /tmp/invocation_command.txt

COPY_INSTRUCTIONS

# Save to tmp with better name
cp "$TEMP_INVOCATION" /tmp/invocation_command.txt
echo "✅ Invocation command saved to: /tmp/invocation_command.txt"
echo ""

# Show how to copy for each platform
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "📋 Copy Command for Your Platform:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

OS=$(uname -s)
case "$OS" in
  Darwin)
    echo "🍎 macOS detected. Copy to clipboard with:"
    echo ""
    echo "  cat /tmp/invocation_command.txt | pbcopy"
    echo ""
    ;;
  Linux)
    echo "🐧 Linux detected. Copy to clipboard with:"
    echo ""
    if command -v xclip &> /dev/null; then
      echo "  cat /tmp/invocation_command.txt | xclip -selection clipboard"
    elif command -v xsel &> /dev/null; then
      echo "  cat /tmp/invocation_command.txt | xsel --clipboard"
    else
      echo "  cat /tmp/invocation_command.txt  # Then manually select and copy"
    fi
    echo ""
    ;;
  *)
    echo "❓ Unknown OS. View and copy manually:"
    echo ""
    echo "  cat /tmp/invocation_command.txt"
    echo ""
    ;;
esac

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "⚡ Quick Start (Once Claude is open):"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "1. In Claude terminal, paste the complete invocation"
echo "2. Press ENTER"
echo "3. Wait 10-15 minutes for Claude to generate all files"
echo "4. Claude will output them as code blocks"
echo "5. Extract and save each file locally"
echo "6. Commit to git: git add -A && git commit -m 'Agent-generated pipeline'"
echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "📚 For detailed instructions, see:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  cat master-invocation-guide.md"
echo ""

rm -f "$TEMP_INVOCATION"

echo ""
echo "✅ Ready! Copy the invocation command and paste into Claude 🚀"
echo ""
