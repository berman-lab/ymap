<?php
	session_start();
	require_once '../constants.php';
	error_reporting(E_ALL | E_STRICT);

	// Only information needed by this script, sent by "js/ajaxfileupload.js" from form defined in "uploader.1.php".
	// This safe data is used to construct the target file location.
	$user             = $_POST["target_user"];
	$genome           = $_POST["target_genome"];
	$project          = $_POST["target_project"];
	if ($genome != "") {
		$target_dir = "../../users/".$user."/genomes/".$genome."/";
	} else if ($project != "") {
		$target_dir = "../../users/".$user."/projects/".$project."/";
	}

	//============================================================================
	// HTML5Uploader ::  Adam Filkor : http://filkor.org
	// Licensed under the MIT license : http://www.opensource.org/licenses/MIT
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	// non-MySQL version.
	//----------------------------------------------------------------------------
	require('UploadHandler.php');
	$upload_handler = new UploadHandler($target_dir);
?>
