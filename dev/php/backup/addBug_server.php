<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	ini_set('display_errors', 1);

	$submitter   = filter_input(INPUT_POST, "submitter",   FILTER_SANITIZE_STRING);
	$description = filter_input(INPUT_POST, "description", FILTER_SANITIZE_STRING);

	if($submitter != "" && $description != ""){
		// Order of replacement
		$order   = array("\r\n", "\n", "\r");
		$replace = '<br>';

		// Processes \r\n's first so they aren't converted twice.
		$description = str_replace($order, $replace, $description);

		$bugtrackFile = $GLOBALS['directory'].$GLOBALS['bugtrackDir'].$GLOBALS['bugtrackFile'];

		// Write the new bug into the tracking file.
		$bugFileHandle = fopen($bugtrackFile, 'a');
		fwrite($bugFileHandle, $submitter."|".$description."|Open|Notes...\n");
		fclose($bugFileHandle);
		chmod($GLOBALS['directory'].$GLOBALS['bugtrackDir'].$GLOBALS['bugtrackFile'], 0777);
	}
?>
