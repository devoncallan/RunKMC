#!/bin/bash

# RunKMC Profiling Script
# Usage: ./profile.sh <input_file> <output_dir> [additional_flags]
# Example: ./profile.sh docs/examples/CRP1_Example.txt output/ --report-polymers --report-sequences

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROFILE_DIR="${SCRIPT_DIR}/profiling"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
PROFILE_LOG="${PROFILE_DIR}/profile_${TIMESTAMP}.txt"
SUMMARY_LOG="${PROFILE_DIR}/summary_${TIMESTAMP}.txt"

# Get RunKMC binary path
BINARY_PATH=$(python3 -c "from runkmc import PATHS; print(PATHS.EXECUTABLE_PATH)" 2>/dev/null)

if [ -z "$BINARY_PATH" ] || [ ! -f "$BINARY_PATH" ]; then
    echo -e "${RED}Error: RunKMC binary not found${NC}"
    echo "Run 'python3 -m runkmc.kmc.build' to build it first"
    exit 1
fi

# Check arguments
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input_file> <output_dir> [additional_flags]"
    echo ""
    echo "Examples:"
    echo "  $0 input.txt output/"
    echo "  $0 input.txt output/ --report-polymers"
    echo "  $0 input.txt output/ --report-polymers --report-sequences"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_DIR="$2"
shift 2
ADDITIONAL_FLAGS="$@"

# Validate input file
if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${RED}Error: Input file not found: ${INPUT_FILE}${NC}"
    exit 1
fi

# Create profiling directory
mkdir -p "${PROFILE_DIR}"

# Clean output directory
rm -rf "${OUTPUT_DIR}"

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  RunKMC Profiling Session${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
echo -e "${GREEN}Configuration:${NC}"
echo "  Binary:       ${BINARY_PATH}"
echo "  Input:        ${INPUT_FILE}"
echo "  Output:       ${OUTPUT_DIR}"
echo "  Flags:        ${ADDITIONAL_FLAGS:-none}"
echo "  Profile log:  ${PROFILE_LOG}"
echo "  Summary log:  ${SUMMARY_LOG}"
echo ""

# Run profiling
echo -e "${YELLOW}Starting profiling (10 second sample)...${NC}"
echo ""

# Start sampling in background
# sample RunKMC 10 -file "${PROFILE_LOG}" &
# SAMPLE_PID=$!
sample RunKMC -file "${PROFILE_LOG}" &  # No time limit
SAMPLE_PID=$!

# Run the simulation
"${BINARY_PATH}" "${INPUT_FILE}" "${OUTPUT_DIR}" ${ADDITIONAL_FLAGS}

# Wait for sampling to complete
wait $SAMPLE_PID 2>/dev/null || true

echo ""
echo -e "${GREEN}✓ Profiling complete${NC}"
echo ""

# Generate summary
echo -e "${YELLOW}Generating summary...${NC}"

python3 << EOF > "${SUMMARY_LOG}"
import re
from collections import defaultdict
from datetime import datetime

print("=" * 80)
print("RUNKMC PROFILING SUMMARY")
print("=" * 80)
print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Input:     ${INPUT_FILE}")
print(f"Output:    ${OUTPUT_DIR}")
print(f"Flags:     ${ADDITIONAL_FLAGS:-none}")
print("=" * 80)
print()

with open('${PROFILE_LOG}', 'r') as f:
    content = f.read()

# Extract function samples
pattern = r'(\d+)\s+([\w:]+(?:\(.*?\))?)\s+\(in'
matches = re.findall(pattern, content)

# Aggregate by function
function_samples = defaultdict(int)
for count, func in matches:
    func = func.split('(')[0].strip()
    function_samples[func] += int(count)

# Sort by sample count
sorted_funcs = sorted(function_samples.items(), key=lambda x: x[1], reverse=True)

print("TOP 20 FUNCTIONS BY CPU TIME")
print("-" * 80)
print(f"{'Samples':<10} {'%':<8} {'Function'}")
print("-" * 80)

total_samples = sum(function_samples.values())
for func, count in sorted_funcs[:20]:
    pct = 100.0 * count / total_samples
    print(f"{count:<10} {pct:>6.1f}%  {func}")

print()
print("=" * 80)
print("CATEGORY BREAKDOWN (Leaf Functions Only)")
print("=" * 80)
print("Note: Excludes parent functions like KMC::run, KMC::step, main, start")
print()

# Exclude parent/wrapper functions that just call other things
exclude_funcs = {'start', 'main', 'KMC::run', 'KMC::step', 'KMC::runToTime',
                 'KMC::updateSystemState', 'SpeciesSet::analyze'}

# Get leaf functions only
leaf_funcs = [(f, c) for f, c in sorted_funcs if f not in exclude_funcs]
leaf_total = sum(c for _, c in leaf_funcs)

# Calculate categories from leaf functions
analysis_time = sum(count for func, count in leaf_funcs if 'analysis' in func.lower() or 'SequenceStats' in func)
malloc_time = sum(count for func, count in leaf_funcs if 'malloc' in func.lower() or 'free' in func.lower() or 'new' in func.lower() or 'szone' in func.lower())
math_time = sum(count for func, count in leaf_funcs if func in ['log', 'exp', 'sqrt', 'pow'])
reaction_time = sum(count for func, count in leaf_funcs if '::react' in func or 'Reaction' in func)
polymer_ops = sum(count for func, count in leaf_funcs if 'Polymer' in func and 'Container' not in func)
container_ops = sum(count for func, count in leaf_funcs if 'Container' in func)
registry_time = sum(count for func, count in leaf_funcs if 'Registry' in func or 'registry' in func)

categories = [
    ("Math (log, exp, etc.)", math_time),
    ("Polymer container operations", container_ops),
    ("Analysis (sequence stats)", analysis_time),
    ("Reaction execution", reaction_time),
    ("Polymer operations", polymer_ops),
    ("Memory allocation/deallocation", malloc_time),
    ("Registry lookups", registry_time),
]

for name, time in sorted(categories, key=lambda x: x[1], reverse=True):
    if time > 0:
        pct = 100.0 * time / total_samples
        bar_length = int(pct / 2)  # Scale to 50 chars max
        bar = "█" * bar_length
        print(f"{name:35} {pct:>5.1f}%  {bar}")

# Show unaccounted time
accounted = sum(t for _, t in categories)
unaccounted = leaf_total - accounted
unaccounted_pct = 100.0 * unaccounted / total_samples
print(f"{'Other/distributed':35} {unaccounted_pct:>5.1f}%  {'█' * int(unaccounted_pct / 2)}")

print()
print("=" * 80)
print("KEY FINDINGS")
print("=" * 80)

# Get top individual functions (excluding parents)
top_leaf = [f for f in leaf_funcs[:5]]

print("Top 5 hottest functions:")
for func, count in top_leaf:
    pct = 100.0 * count / total_samples
    print(f"  • {func:45} {pct:>5.1f}%")

print()

# Identify top category bottleneck
top_category = max(categories, key=lambda x: x[1])
pct = 100.0 * top_category[1] / total_samples

if pct > 10:
    print(f"🔥 PRIMARY CATEGORY: {top_category[0]} ({pct:.1f}% of runtime)")
    print()

    if 'math' in top_category[0].lower():
        print("   Analysis:")
        print("   → Mathematical operations (log, exp) are expected in KMC")
        print("   → Gillespie algorithm requires log() per step")
        print("   → This is likely optimal unless you switch algorithms")
    elif 'container' in top_category[0].lower():
        print("   Optimization opportunities:")
        print("   → Profile removeRandomPolymer() - may need better data structure")
        print("   → Consider using std::unordered_set or vector with swap-and-pop")
        print("   → Benchmark different random selection strategies")
    elif 'analysis' in top_category[0].lower():
        print("   Optimization opportunities:")
        print("   → Sequence analysis has O(n) complexity per polymer")
        print("   → Consider making --report-sequences less frequent")
        print("   → Pre-allocate SequenceStats objects")
    elif 'malloc' in top_category[0].lower():
        print("   Optimization opportunities:")
        print("   → Excessive memory allocation detected")
        print("   → Consider object pooling for Polymer objects")
        print("   → Pre-allocate containers with reserve()")
    elif 'polymer' in top_category[0].lower():
        print("   Optimization opportunities:")
        print("   → Profile Polymer::terminate() and insertPolymer()")
        print("   → Consider more efficient polymer data structures")
else:
    print("✓ No single dominant bottleneck found")
    print("  Runtime is well-distributed across components")
    print()
    print("  This is good! It means:")
    print("  • No obvious low-hanging fruit for optimization")
    print("  • The simulation is well-balanced")
    print("  • Any speedup would require algorithmic changes")

print()
print("Full profile saved to: ${PROFILE_LOG}")
print("=" * 80)
EOF

# Display summary to terminal
cat "${SUMMARY_LOG}"

echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  Profiling Results Saved${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "  Full profile: ${PROFILE_LOG}"
echo "  Summary:      ${SUMMARY_LOG}"
echo ""
echo "To view full profile:"
echo "  open ${PROFILE_LOG}"
echo ""
echo "To list all profiles:"
echo "  ls -lht ${PROFILE_DIR}/"
echo ""
