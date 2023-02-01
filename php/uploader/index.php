<?php
	session_start();
	error_reporting(E_ALL);  // | E_STRICT
	require_once '../../constants.php';
	require_once '../../POST_validation.php';
	ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.
	if(!isset($_SESSION['logged_on'])){
		session_destroy();
		header('Location: ../../');
	}

	// Load user string from session.
	$user    = $_SESSION['user'];

	// Sanitize input strings.
	$genome  = sanitize_POST("target_genome");
	$project = sanitize_POST("target_project");

	// Either a genome or project string should be present as something other than "".
	if ($genome != "") {
		// Confirm if requested genome exists.
		$genome_dir = "../../users/".$user."/genomes/".$genome;
		if (!is_dir($genome_dir)) {
			// Genome doesn't exist, should never happen: Force logout.
			session_destroy();
			header('Location: ../../');
		} else {
			$target_dir = $genome_dir."/";
		}
	} else if ($project != "") {
		// Confirm if requested project exists.
		$project_dir = "../../users/".$user."/projects/".$project;
		if (!is_dir($project_dir)) {
			// Project doesn't exist, should never happen: Force logout.
			session_destroy();
			header('Location: ../../');
		} else {
			$target_dir = $project_dir."/";
		}
	}

	// Data received from: "js/ajaxfileupload.js" from form defined in "uploader.1.php".


	//============================================================================
	// HTML5Uploader ::  Adam Filkor : http://filkor.org
	// Licensed under the MIT license : http://www.opensource.org/licenses/MIT
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	// non-MySQL version.
	//----------------------------------------------------------------------------
	require('UploadHandler.php');
	$upload_handler = new UploadHandler($target_dir);
?>
