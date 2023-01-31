<?php
	session_start();
	error_reporting(E_ALL);  // | E_STRICT
	require_once '../constants.php';
	ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.
	if(!isset($_SESSION['logged_on'])){
		session_destroy();
		header('Location: user.login.php');
	}

	// Load user string from session.
	$user    = $_SESSION['user'];

	// Sanitize input strings.
	$genome  = trim(filter_input(INPUT_POST, "target_genome", FILTER_SANITIZE_STRING));	// strip out any html tags.
	$genome  = str_replace(" ","_",$genome);						// convert any spaces to underlines.
	$genome  = preg_replace("/[\s\W]+/", "", $genome);					// remove everything but alphanumeric characters and underlines.
	$project = trim(filter_input(INPUT_POST, "target_project", FILTER_SANITIZE_STRING));
	$project = str_replace(" ","_",$project);
	$project = preg_replace("/[\s\W]+/", "", $project);


	// Either a genome or project string should be present as something other than "".
	if ($genome != "") {
		// Confirm if requested genome exists.
		$genome_dir = "../../users/".$user."/genomes/".$genome;
		if (!is_dir($genome_dir)) {
			// Genome doesn't exist, should never happen: Force logout.
			session_destroy();
			header('Location: user.login.php');
		} else {
			$target_dir = $genome_dir."/";
		}
	} else if ($project != "") {
		// Confirm if requested project exists.
		$project_dir = "../../users/".$user."/projects/".$project;
		if (!is_dir($project_dir)) {
			// Project doesn't exist, should never happen: Force logout.
			session_destroy();
			header('Location: user.login.php');
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
