<?php
	session_start();
	error_reporting(E_ALL | E_STRICT);

	// Only information needed by this script, sent by "js/ajaxfileupload.js" from form defined in "uploader.1.php".
	// This safe data is used to construct the target file location.

	$bad_chars = array("~","@","#","$","%","^","&","*","(",")","+","=","|","{","}","<",">","?",".",",","\\","/","'",'"',"[","]","!");
	$user     = str_replace($bad_chars,"",trim(filter_input(INPUT_POST, "target_user",    FILTER_SANITIZE_STRING)));
	$genome   = str_replace($bad_chars,"",trim(filter_input(INPUT_POST, "target_genome",  FILTER_SANITIZE_STRING)));
	$project  = str_replace($bad_chars,"",trim(filter_input(INPUT_POST, "target_project", FILTER_SANITIZE_STRING)));

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
