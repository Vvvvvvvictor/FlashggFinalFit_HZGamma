#!/bin/bash

# Script to automatically optimize background fit range (mgg_low and mgg_high)
# Using simplified binary search algorithm to find the optimal range
# Low range: 95-110, High range: 150-180
# Author: Jie Han (Modified with simplified binary search algorithm)

source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
source setup.sh
cd Background
# make clean; make -j 16

# Set default values
DATA_FILE="/eos/home-j/jiehan/root/input_finalfit/background/ws/output_Data_all.root"
OUTDIR="outdir_fiducialAnalysis"
ISDATA=1
YEAR="all"
DATAEXT="-Data"
CATEGORIES=("VBF0" "VBF1" "VBF2" "VBF3") # "VBF0" "VBF1" "VBF2" "VBF3"

# Define category-specific parameters
# Significance targets
declare -A SIGNIFICANCE_TARGETS=(
  ["VBF0"]=0.286
  ["VBF1"]=0.588
  ["VBF2"]=0.939
  ["VBF3"]=0.948
)

# Low range boundaries for each category
declare -A LOW_MIN_RANGES=(
  ["VBF0"]=100
  ["VBF1"]=97
  ["VBF2"]=91
  ["VBF3"]=91
)
declare -A LOW_MAX_RANGES=(
  ["VBF0"]=110
  ["VBF1"]=105
  ["VBF2"]=100
  ["VBF3"]=100
)

# High range boundaries for each category
declare -A HIGH_MIN_RANGES=(
  ["VBF0"]=160
  ["VBF1"]=160
  ["VBF2"]=150
  ["VBF3"]=150
)
declare -A HIGH_MAX_RANGES=(
  ["VBF0"]=180
  ["VBF1"]=180
  ["VBF2"]=165
  ["VBF3"]=165
)
DEFAULT_SIGNIFICANCE_TARGET=1.0
SIGNIFICANCE_TOLERANCE=0.15  # Tolerance for significance (15%)

# Output directory for results
RESULTS_DIR="${OUTDIR}/range_optimization"
mkdir -p "${RESULTS_DIR}"

# Output summary file
SUMMARY_FILE="${RESULTS_DIR}/range_optimization_summary.json"
LOG_FILE="${RESULTS_DIR}/binary_search_progress.log"

# Initialize JSON output file
echo '{' > "${SUMMARY_FILE}"
echo '  "range_results": {' >> "${SUMMARY_FILE}"

# Function to extract significance value from combine output
extract_significance() {
  local logfile=$1
  if [ -f "$logfile" ]; then
    sig=$(grep "Significance:" "$logfile" | awk '{print $2}')
    if [ -z "$sig" ]; then
      echo "0.0"
    else
      echo $sig
    fi
  else
    echo "0.0"
  fi
}

# Function to evaluate a specific range configuration
evaluate_range() {
  local cat=$1
  local cat_idx=$2
  local mgg_low=$3
  local mgg_high=$4
  local range_dir=$5
  
  # Ensure integers
  mgg_low=$((mgg_low))
  mgg_high=$((mgg_high))
  
  echo "Testing range: ${mgg_low}-${mgg_high}" >> "${LOG_FILE}"
  
  # Create output directory for this range
  mkdir -p "${range_dir}"
  
  # Run F-test with this range
  echo "Running F-test..." >> "${LOG_FILE}"
  echo "F-test command: ./bin/fTest -i $DATA_FILE --saveMultiPdf $OUTDIR/CMS-HGG_multipdf_${cat}_${mgg_low}_${mgg_high}.root -D ${range_dir}/bkgfTest_${cat}_${mgg_low}_${mgg_high} -f ${cat} --mgg_low ${mgg_low} --mgg_high ${mgg_high} --isData ${ISDATA} --year ${YEAR} --catOffset ${cat_idx} --blindFit" >> "${LOG_FILE}"
  ./bin/fTest -i $DATA_FILE --saveMultiPdf $OUTDIR/CMS-HGG_multipdf_${cat}_${mgg_low}_${mgg_high}.root \
    -D ${range_dir}/bkgfTest_${cat}_${mgg_low}_${mgg_high} \
    -f ${cat} --mgg_low ${mgg_low} --mgg_high ${mgg_high} \
    --isData ${ISDATA} --year ${YEAR} --catOffset ${cat_idx} --blindFit > /dev/null 2>&1

  echo "Finished F-test for range ${mgg_low}-${mgg_high}" >> "${LOG_FILE}"
  
  # Check if chi2 results file exists
  CHI2_FILE="${range_dir}/bkgfTest_${cat}_${mgg_low}_${mgg_high}/fTest_chi2_${cat}.json"
  if [ -f "${CHI2_FILE}" ]; then
    # Extract average chi2 and function count from JSON
    chi2_avg=$(grep -o '"average_chi2_per_ndof":[^,}]*' "${CHI2_FILE}" | cut -d: -f2)
    func_count=$(grep -o '"num_functions":[^,}]*' "${CHI2_FILE}" | cut -d: -f2)
    
    # Handle empty values
    # Check if chi2_avg is nan
    if [[ "$chi2_avg" == *"nan"* ]]; then
      chi2_avg="999"
    fi
    if [ -z "$chi2_avg" ]; then chi2_avg="999"; fi
    if [ -z "$func_count" ]; then func_count="0"; fi
    
    echo "Average chi2/ndof: ${chi2_avg}, Function count: ${func_count}" >> "${LOG_FILE}"
    
    # Generate workspace and run significance
    echo "Running Text2Workspace and significance calculation..." >> "${LOG_FILE}"
    cd ../Combine
    cp ../Datacard/Datacard_${cat}.txt .
    # rm -rf Models 2>/dev/null
    mkdir -p Models/background
    cp ../Background/${OUTDIR}/CMS-HGG_multipdf_${cat}_${mgg_low}_${mgg_high}.root Models/background/CMS-HGG_multipdf_${cat}.root
    mkdir -p Models/signal
    cp -r ../Signal/outdir_packaged/*.root Models/signal/
    
    # Run Text2Workspace and calculate significance
    SIGNIFICANCE_LOG="../Background/${range_dir}/significance.log"
    python RunText2Workspace.py --ext _${cat} --mode mu_fiducial --batch local
    combine Datacard_${cat}_mu_fiducial.root -M Significance -t -1 --expectSignal=1 -m 125.0 -n _${cat} > "${SIGNIFICANCE_LOG}" 2>&1
    significance=$(extract_significance "${SIGNIFICANCE_LOG}")
    
    # Handle empty significance
    if [ -z "$significance" ]; then significance="0.0"; fi
    
    cd ../Background
    echo "Significance: ${significance}" >> "${LOG_FILE}"
    
    # Only output the necessary values (without additional text)
    # 添加标记1表示成功生成chi2结果
    echo "${chi2_avg} ${func_count} ${significance} 1"
  else
    echo "ERROR: Chi2 results file not found for range ${mgg_low}-${mgg_high}" >> "${LOG_FILE}"
    # 添加标记0表示未成功生成chi2结果
    echo "0 0 0 0"
  fi
}

# Function to record results to JSON file
record_results() {
  local cat=$1
  local mgg_low=$2
  local mgg_high=$3
  local chi2_avg=$4
  local func_count=$5
  local significance=$6
  local is_comma=$7
  
  # Ensure integers
  mgg_low=$((mgg_low))
  mgg_high=$((mgg_high))
  
  # Handle empty values with default 0
  if [ -z "$chi2_avg" ]; then chi2_avg="0"; fi
  if [ -z "$func_count" ]; then func_count="0"; fi
  if [ -z "$significance" ]; then significance="0"; fi
  
  # Add comma if needed
  if [ "$is_comma" = "true" ]; then
    echo ',' >> "${SUMMARY_FILE}"
  fi
  
  # Write properly formatted JSON
  echo "    \"${mgg_low}_${mgg_high}\": {" >> "${SUMMARY_FILE}"
  echo "      \"chi2_avg\": ${chi2_avg}," >> "${SUMMARY_FILE}"
  echo "      \"func_count\": ${func_count}," >> "${SUMMARY_FILE}"
  echo "      \"significance\": ${significance}" >> "${SUMMARY_FILE}"
  echo "    }" >> "${SUMMARY_FILE}"
}

# Function to compare two range configurations and determine which is better
# Returns 0 if first range is better, 1 if second range is better
compare_ranges() {
  local chi2_1=$1
  local func_1=$2
  local sig_1=$3
  local success_1=$4
  local chi2_2=$5
  local func_2=$6
  local sig_2=$7
  local success_2=$8
  local sig_target=$9
  
  # 首先检查是否有成功生成chi2结果，如果一个成功而另一个失败，优先选择成功的
  if [ "$success_1" -eq 1 ] && [ "$success_2" -eq 0 ]; then
    # 第一个成功生成了chi2结果，而第二个没有，选择第一个
    echo "选择${chi2_1} ${func_1} ${sig_1}，因为它成功生成了chi2结果而另一个没有" >> "${LOG_FILE}"
    return 0
  elif [ "$success_1" -eq 0 ] && [ "$success_2" -eq 1 ]; then
    # 第二个成功生成了chi2结果，而第一个没有，选择第二个
    echo "选择${chi2_2} ${func_2} ${sig_2}，因为它成功生成了chi2结果而另一个没有" >> "${LOG_FILE}"
    return 1
  elif [ "$success_1" -eq 0 ] && [ "$success_2" -eq 0 ]; then
    # 两个都没有成功生成chi2结果，随便返回一个（这种情况应该避免）
    echo "警告：两个配置都未能生成chi2结果，随机选择" >> "${LOG_FILE}"
    return 0
  fi
  
  # 以下是两个都成功生成chi2结果的情况
  # Handle empty values with default 0
  if [ -z "$chi2_1" ]; then chi2_1="100"; fi
  if [ -z "$func_1" ]; then func_1="0"; fi
  if [ -z "$sig_1" ]; then sig_1="0"; fi
  if [ -z "$chi2_2" ]; then chi2_2="100"; fi
  if [ -z "$func_2" ]; then func_2="0"; fi
  if [ -z "$sig_2" ]; then sig_2="0"; fi
  
  # Calculate absolute difference from target significance
  local sig_diff_1=$(echo "${sig_1} - ${sig_target}" | bc -l)
  local sig_diff_1_abs=$(echo "${sig_diff_1#-}" | bc -l) # absolute value
  local sig_diff_2=$(echo "${sig_2} - ${sig_target}" | bc -l)
  local sig_diff_2_abs=$(echo "${sig_diff_2#-}" | bc -l) # absolute value
  
  local within_tolerance_1=$(echo "${sig_diff_1_abs} <= ${SIGNIFICANCE_TOLERANCE}" | bc -l)
  local within_tolerance_2=$(echo "${sig_diff_2_abs} <= ${SIGNIFICANCE_TOLERANCE}" | bc -l)
  
  # If both within tolerance, compare chi2 and function count
  if [ "$within_tolerance_1" -eq 1 ] && [ "$within_tolerance_2" -eq 1 ]; then
    # Prefer lower chi2, then higher function count
    if [ ${func_1} -gt ${func_2} ]; then
      return 0
    elif [ $(echo "${chi2_1} < ${chi2_2}" | bc -l) -eq 1 ]; then
      return 0
    elif [ $(echo "${chi2_1} > ${chi2_2}" | bc -l) -eq 1 ]; then
      return 1
    else
      return 1
    fi
  # If only one is within tolerance, prefer that one
  elif [ "$within_tolerance_1" -eq 1 ]; then
    return 0
  elif [ "$within_tolerance_2" -eq 1 ]; then
    return 1
  # If neither is within tolerance, prefer higher significance
  elif [ $(echo "${sig_1} > ${sig_2}" | bc -l) -eq 1 ]; then
    return 0
  else
    return 1
  fi
}

# Default values (used if category not found in arrays)
DEFAULT_LOW_MIN=95
DEFAULT_LOW_MAX=110
DEFAULT_HIGH_MIN=150
DEFAULT_HIGH_MAX=180

# Initialize log file
echo "Binary Search Optimization Log - $(date)" > "${LOG_FILE}"

# Process each category
first_cat=true
for cat_idx in ${!CATEGORIES[@]}; do
  cat=${CATEGORIES[$cat_idx]}
  # Get the significance target for this category
  SIGNIFICANCE_TARGET=${SIGNIFICANCE_TARGETS[$cat]:-$DEFAULT_SIGNIFICANCE_TARGET}
  
  echo "=========================================" >> "${LOG_FILE}"
  echo "Processing category: ${cat}" | tee -a "${LOG_FILE}"
  echo "Significance target: ${SIGNIFICANCE_TARGET}" >> "${LOG_FILE}"
  echo "=========================================" >> "${LOG_FILE}"
  
  # Start the category entry in JSON
  if [ "$first_cat" = true ]; then
    first_cat=false
  else
    echo "," >> "${SUMMARY_FILE}"
  fi
  
  echo "    \"${cat}\": {" >> "${SUMMARY_FILE}"
  
  # Get range boundaries for this category
  LOW_MIN=${LOW_MIN_RANGES[$cat]:-$DEFAULT_LOW_MIN}
  LOW_MAX=${LOW_MAX_RANGES[$cat]:-$DEFAULT_LOW_MAX}
  HIGH_MIN=${HIGH_MIN_RANGES[$cat]:-$DEFAULT_HIGH_MIN}
  HIGH_MAX=${HIGH_MAX_RANGES[$cat]:-$DEFAULT_HIGH_MAX}
  
  echo "Search ranges for ${cat} - Low: ${LOW_MIN}-${LOW_MAX}, High: ${HIGH_MIN}-${HIGH_MAX}" >> "${LOG_FILE}"
  
  # Best configuration tracking
  best_chi2_avg=999
  best_func_count=0
  best_significance=0
  best_low=$LOW_MAX
  best_high=$HIGH_MIN
  
  # First configuration flag for JSON
  first_config=true
  
  # Phase 1: Optimize the low boundary with fixed high boundary
  echo "Phase 1: Optimizing low boundary..." | tee -a "${LOG_FILE}"
  
  low_min=$LOW_MIN
  low_max=$LOW_MAX
  high_fixed=$((HIGH_MIN + (HIGH_MAX - HIGH_MIN) / 5))
  
  # First evaluate low and high points
  RANGE_DIR="${RESULTS_DIR}/${cat}_${low_min}_${high_fixed}"
  low_results=($(evaluate_range $cat $cat_idx $low_min $high_fixed "$RANGE_DIR"))
  while [ ${low_results[3]} -eq 0 ]; do
    echo "Retrying low range evaluation..." >> "${LOG_FILE}"
    low_min=$((low_min + 1))
    RANGE_DIR="${RESULTS_DIR}/${cat}_${low_min}_${high_fixed}"
    low_results=($(evaluate_range $cat $cat_idx $low_min $high_fixed "$RANGE_DIR"))
  done
  
  record_results $cat $low_min $high_fixed ${low_results[0]} ${low_results[1]} ${low_results[2]} $first_config
  first_config=false
  
  RANGE_DIR="${RESULTS_DIR}/${cat}_${low_max}_${high_fixed}"
  high_results=($(evaluate_range $cat $cat_idx $low_max $high_fixed "$RANGE_DIR"))
  while [ ${high_results[3]} -eq 0 ]; do
    echo "Retrying high range evaluation..." >> "${LOG_FILE}"
    low_max=$((low_max - 1))
    RANGE_DIR="${RESULTS_DIR}/${cat}_${low_max}_${high_fixed}"
    high_results=($(evaluate_range $cat $cat_idx $low_max $high_fixed "$RANGE_DIR"))
  done
  
  record_results $cat $low_max $high_fixed ${high_results[0]} ${high_results[1]} ${high_results[2]} true
  
  # Main binary search loop for low boundary
  while [ $((low_max - low_min)) -gt 1 ]; do
    echo "Low search range: ${low_min}-${low_max}" | tee -a "${LOG_FILE}"
    
    # Calculate midpoint
    low_mid=$(( (low_min + low_max) / 2 ))

    RANGE_DIR="${RESULTS_DIR}/${cat}_${low_mid}_${high_fixed}"
    mid_results=($(evaluate_range $cat $cat_idx $low_mid $high_fixed "$RANGE_DIR"))

    while [ ${mid_results[3]} -eq 0 ] && [$(low_mid -gt low_min)] && [$(low_max -gt low_mid)] ; do
      if [ $high_better_than_low -eq 0 ]; then
        low_mid=$((low_mid + 1))
      else
        low_mid=$((low_mid - 1))
      fi
      mid_results=($(evaluate_range $cat $cat_idx $low_mid $high_fixed "$RANGE_DIR"))
    done
    
    record_results $cat $low_mid $high_fixed ${mid_results[0]} ${mid_results[1]} ${mid_results[2]} true
    
    # Compare and update search space based on which is better
    compare_ranges ${high_results[0]} ${high_results[1]} ${high_results[2]} ${high_results[3]} \
                  ${low_results[0]} ${low_results[1]} ${low_results[2]} ${low_results[3]} \
                  $SIGNIFICANCE_TARGET
    high_better_than_low=$?
    
    if [ $high_better_than_low -eq 0 ]; then
      # If high is better than low, replace low with mid
      low_min=$low_mid
      low_results=("${mid_results[@]}")
    else
      # If low is better than high, replace high with mid
      low_max=$low_mid
      high_results=("${mid_results[@]}")
    fi
  done
  
  # Determine best low value
  compare_ranges ${low_results[0]} ${low_results[1]} ${low_results[2]} ${low_results[3]} \
               ${high_results[0]} ${high_results[1]} ${high_results[2]} ${high_results[3]} \
               $SIGNIFICANCE_TARGET
  if [ $? -eq 0 ]; then
    best_low=$low_min
    best_chi2_avg=${low_results[0]}
    best_func_count=${low_results[1]}
    best_significance=${low_results[2]}
  else
    best_low=$low_max
    best_chi2_avg=${high_results[0]}
    best_func_count=${high_results[1]}
    best_significance=${high_results[2]}
  fi
  
  echo "Best low value found: ${best_low}" | tee -a "${LOG_FILE}"
  
  # Phase 2: Optimize the high boundary with fixed optimized low boundary
  echo "Phase 2: Optimizing high boundary..." | tee -a "${LOG_FILE}"
  
  high_min=$HIGH_MIN
  high_max=$HIGH_MAX
  low_fixed=$best_low
  
  # First evaluate low and high points
  RANGE_DIR="${RESULTS_DIR}/${cat}_${low_fixed}_${high_min}"
  low_results=($(evaluate_range $cat $cat_idx $low_fixed $high_min "$RANGE_DIR"))
  
  record_results $cat $low_fixed $high_min ${low_results[0]} ${low_results[1]} ${low_results[2]} true
  
  RANGE_DIR="${RESULTS_DIR}/${cat}_${low_fixed}_${high_max}"
  high_results=($(evaluate_range $cat $cat_idx $low_fixed $high_max "$RANGE_DIR"))
  
  record_results $cat $low_fixed $high_max ${high_results[0]} ${high_results[1]} ${high_results[2]} true
  
  # Main binary search loop for high boundary
  while [ $((high_max - high_min)) -gt 1 ]; do
    echo "High search range: ${high_min}-${high_max}" | tee -a "${LOG_FILE}"
    
    # Compare and update search space based on which is better
    compare_ranges ${high_results[0]} ${high_results[1]} ${high_results[2]} ${high_results[3]} \
                  ${low_results[0]} ${low_results[1]} ${low_results[2]} ${low_results[3]} \
                  $SIGNIFICANCE_TARGET
    high_better_than_low=$?

    # Calculate midpoint
    high_mid=$(( (high_min + high_max) / 2 ))

    # Evaluate midpoint
    RANGE_DIR="${RESULTS_DIR}/${cat}_${low_fixed}_${high_mid}"
    mid_results=($(evaluate_range $cat $cat_idx $low_fixed $high_mid "$RANGE_DIR"))
    while [ ${mid_results[3]} -eq 0] && [$((high_max - high_mid)) -gt 1] && [$((high_mid - high_min)) -gt 1 ]; do
      if [ $high_better_than_low -eq 0 ]; then
        high_mid=$((high_mid + 1))
      else
        high_mid=$((high_mid - 1))
      fi
    done    

    record_results $cat $low_fixed $high_mid ${mid_results[0]} ${mid_results[1]} ${mid_results[2]} true
    
    if [ $high_better_than_low -eq 0 ]; then
      # If high is better than low, replace low with mid
      high_min=$high_mid
      low_results=("${mid_results[@]}")
    else
      # If low is better than high, replace high with mid
      high_max=$high_mid
      high_results=("${mid_results[@]}")
    fi
  done
  
  # Determine best high value
  compare_ranges ${low_results[0]} ${low_results[1]} ${low_results[2]} ${low_results[3]} \
               ${high_results[0]} ${high_results[1]} ${high_results[2]} ${high_results[3]} \
               $SIGNIFICANCE_TARGET
  if [ $? -eq 0 ]; then
    best_high=$high_min
    best_chi2_avg=${low_results[0]}
    best_func_count=${low_results[1]}
    best_significance=${low_results[2]}
  else
    best_high=$high_max
    best_chi2_avg=${high_results[0]}
    best_func_count=${high_results[1]}
    best_significance=${high_results[2]}
  fi
  
  echo "Best high value found: ${best_high}" | tee -a "${LOG_FILE}"
  
  # Add best configuration summary
  echo "      \"best_config\": {" >> "${SUMMARY_FILE}"
  echo "        \"mgg_low\": ${best_low}," >> "${SUMMARY_FILE}"
  echo "        \"mgg_high\": ${best_high}," >> "${SUMMARY_FILE}" 
  echo "        \"chi2_avg\": ${best_chi2_avg}," >> "${SUMMARY_FILE}"
  echo "        \"func_count\": ${best_func_count}," >> "${SUMMARY_FILE}"
  echo "        \"significance\": ${best_significance}" >> "${SUMMARY_FILE}"
  echo "      }" >> "${SUMMARY_FILE}"
  
  # Close the category entry
  echo "    }" >> "${SUMMARY_FILE}"
  
  echo "============= Final Results =============" | tee -a "${LOG_FILE}"
  echo "Best configuration for ${cat}: range=${best_low}-${best_high}" | tee -a "${LOG_FILE}"
  echo "Significance: ${best_significance}, chi2_avg: ${best_chi2_avg}, func_count: ${best_func_count}" | tee -a "${LOG_FILE}"
  echo "=========================================" | tee -a "${LOG_FILE}"
  
  # Copy the best model to the main output directory
  cp ${OUTDIR}/CMS-HGG_multipdf_${cat}_${best_low}_${best_high}.root ${OUTDIR}/CMS-HGG_multipdf_${cat}.root
  
  # Save the best range to a JSON file for each category
  BEST_RANGE_FILE="${OUTDIR}/best_range_${cat}.json"
  echo "{" > "${BEST_RANGE_FILE}"
  echo "  \"mgg_low\": ${best_low}," >> "${BEST_RANGE_FILE}"
  echo "  \"mgg_high\": ${best_high}," >> "${BEST_RANGE_FILE}"
  echo "  \"chi2_avg\": ${best_chi2_avg}," >> "${BEST_RANGE_FILE}"
  echo "  \"func_count\": ${best_func_count}," >> "${BEST_RANGE_FILE}"
  echo "  \"significance\": ${best_significance}" >> "${BEST_RANGE_FILE}"
  echo "}" >> "${BEST_RANGE_FILE}"
done

# Close the JSON structure
echo '  },' >> "${SUMMARY_FILE}"

# Add timestamp
echo '  "timestamp": "'"$(date '+%Y-%m-%d %H:%M:%S')"'"' >> "${SUMMARY_FILE}"
echo "}" >> "${SUMMARY_FILE}"

echo "Optimization complete! Results saved to ${SUMMARY_FILE}" | tee -a "${LOG_FILE}"
echo "Best configurations have been copied to the main output directory."

# Generate the final workspace files for all categories
echo "Generating final workspaces for all categories..." | tee -a "${LOG_FILE}"

# Run combine for each category with best range
cd ../Combine
for cat in "${CATEGORIES[@]}"; do
  echo "Running significance test for ${cat} with optimized range..." | tee -a "${LOG_FILE}"
  cp ../Datacard/Datacard_${cat}.txt .
  rm -rf Models 2>/dev/null
  mkdir -p Models/background Models/signal
  cp ../${OUTDIR}/CMS-HGG_multipdf_${cat}.root Models/background/
  cp -r ../Models/signal Models/
  
  python RunText2Workspace.py --ext _${cat} --mode mu_fiducial --batch local
  combine Datacard_${cat}_mu_fiducial.root -M Significance -t -1 --expectSignal=1 -m 125.0 -n _${cat}_final
done

echo "Optimization and workspace generation complete!" | tee -a "${LOG_FILE}"
exit 0
