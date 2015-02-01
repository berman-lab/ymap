<?php
	session_start();
	require_once 'constants.php';
	error_reporting(E_ALL | E_STRICT);

	$target_dir        = $_POST["target_dir"];
	$conclusion_script = $_POST["conclusion_script"];
	$key               = $_POST["key"];

	//============================================================================
	// HTML5Uploader ::  Adam Filkor : http://filkor.org
	// Licensed under the MIT license : http://www.opensource.org/licenses/MIT
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	// non-MySQL version.
	//----------------------------------------------------------------------------
	require('UploadHandler.php');
	$upload_handler = new UploadHandler($target_dir, $conclusion_script, $key);
?>
