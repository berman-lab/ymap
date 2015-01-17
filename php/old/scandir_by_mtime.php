<?php
function scandir_by_mtime($folder) {
	$dircontent = scandir($folder);
	$arr = array();
	foreach($dircontent as $filename) {
		if ($filename != '.' && $filename != '..') {
			if (filemtime($folder.$filename) === false) return false;
			$dat = date("YmdHis", filemtime($folder.$filename));
			$arr[$dat] = $filename;
		}
	}
	if (!ksort($arr)) return false;
	return $arr;
}
?>
