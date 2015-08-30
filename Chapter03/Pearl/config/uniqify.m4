# **************************************************************************
# AC_UNIQIFY_LIST( VARIABLE, LIST )
#
# This macro filters out redundant items from a list.  
#
AC_DEFUN([AC_UNIQIFY_LIST], [
ac_uniqued_list=
for ac_item in $2; do
  eval ac_eval_item="$ac_item"
  eval ac_eval_item="$ac_eval_item"
  if test x"$ac_uniqued_list" = x; then
    ac_uniqued_list="$ac_item"
  else
    ac_unique=true
    for ac_uniq in $ac_uniqued_list; do
      eval ac_eval_uniq="$ac_uniq"
      eval ac_eval_uniq="$ac_eval_uniq"
      test x"$ac_eval_item" = x"$ac_eval_uniq" && ac_unique=false
    done
    $ac_unique && ac_uniqued_list="$ac_uniqued_list $ac_item"
  fi
done
$1=$ac_uniqued_list
# unset ac_eval_item
# unset ac_eval_uniq
]) # AC_UNIQIFY_LIST

