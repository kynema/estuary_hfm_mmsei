#!/bin/bash
#
# cleanup_old_checkpoints.sh
#
# Removes old checkpoint (chk*) and plot (plt*) directories, keeping only the highest numbered one.
# If a directory has only one chk or one plt, nothing is removed.
#
# Usage:
#   ./cleanup_old_checkpoints.sh          # Dry-run: shows what would be deleted
#   ./cleanup_old_checkpoints.sh --delete # Actually delete old files
#

DRY_RUN=true

if [ "$1" = "--delete" ]; then
  DRY_RUN=false
  echo "WARNING: Running in DELETE mode. Old checkpoint and plot directories will be removed."
  read -p "Are you sure? (yes/no): " confirm
  if [ "$confirm" != "yes" ]; then
    echo "Cancelled."
    exit 0
  fi
fi

if [ "$DRY_RUN" = true ]; then
  echo "========================================================================="
  echo "DRY RUN MODE: Showing what would be deleted"
  echo "========================================================================="
else
  echo "========================================================================="
  echo "DELETE MODE: Removing old checkpoint and plot directories"
  echo "========================================================================="
fi

total_dirs=0
total_chk_removed=0
total_plt_removed=0

for dir in flat-drag-* slanted-drag-*; do
  [ -d "$dir" ] || continue
  
  # Get sorted lists of chk and plt directories
  chk_list=$(ls -d "$dir"/chk* 2>/dev/null | sort -V)
  plt_list=$(ls -d "$dir"/plt* 2>/dev/null | sort -V)
  
  # Count lines
  chk_count=$(echo "$chk_list" | grep -c . || true)
  plt_count=$(echo "$plt_list" | grep -c . || true)
  
  # Skip if empty
  if [ -z "$chk_list" ]; then chk_count=0; fi
  if [ -z "$plt_list" ]; then plt_count=0; fi
  
  if [ "$chk_count" -gt 1 ] || [ "$plt_count" -gt 1 ]; then
    total_dirs=$((total_dirs + 1))
    
    # Prepare output string
    output="$dir: deleting"
    remaining_files=""
    
    if [ "$chk_count" -gt 1 ]; then
      # Get the last item (highest number) - this is what remains
      last_chk=$(echo "$chk_list" | tail -1)
      last_chk_name=$(basename "$last_chk")
      remaining_files="$last_chk_name"
      
      output="$output $((chk_count - 1))/$chk_count chk files,"
      
      # Delete old ones
      old_chks=$(echo "$chk_list" | head -n $((chk_count - 1)))
      while IFS= read -r old_dir; do
        [ -n "$old_dir" ] || continue
        if [ "$DRY_RUN" = false ]; then
          rm -rf "$old_dir"
          total_chk_removed=$((total_chk_removed + 1))
        else
          total_chk_removed=$((total_chk_removed + 1))
        fi
      done <<< "$old_chks"
    fi
    
    if [ "$plt_count" -gt 1 ]; then
      # Get the last item (highest number) - this is what remains
      last_plt=$(echo "$plt_list" | tail -1)
      last_plt_name=$(basename "$last_plt")
      [ -n "$remaining_files" ] && remaining_files="$remaining_files, "
      remaining_files="$remaining_files$last_plt_name"
      
      output="$output $((plt_count - 1))/$plt_count plt files."
      
      # Delete old ones
      old_plts=$(echo "$plt_list" | head -n $((plt_count - 1)))
      while IFS= read -r old_dir; do
        [ -n "$old_dir" ] || continue
        if [ "$DRY_RUN" = false ]; then
          rm -rf "$old_dir"
          total_plt_removed=$((total_plt_removed + 1))
        else
          total_plt_removed=$((total_plt_removed + 1))
        fi
      done <<< "$old_plts"
    fi
    
    echo -e "$output\n   Remaining: $remaining_files"
  fi
done

echo ""
echo "========================================================================="
if [ "$DRY_RUN" = true ]; then
  echo "DRY RUN SUMMARY: Found $total_dirs directories with old checkpoints/plots"
  echo "To actually delete these directories, run:"
  echo "  ./cleanup_old_checkpoints.sh --delete"
else
  echo "DELETION COMPLETE"
  echo "Removed $total_chk_removed old CHK directories"
  echo "Removed $total_plt_removed old PLT directories"
fi
echo "========================================================================="""
