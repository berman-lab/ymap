<?php
	session_start();
        error_reporting(E_ALL);
        require_once '../../constants.php';
        require_once '../../POST_validation.php';
        ini_set('display_errors', 1);

        // If the user is not logged on, redirect to login page.
        if(!isset($_SESSION['logged_on'])){
                session_destroy();
                header('Location: .');
        }

	// Load user string from session.
	$user    = $_SESSION['user'];


// Only information needed by this script, sent by "js/ajaxfileupload.js" from form defined in "uploader.1.php".
// This safe data is used to construct the target file locations.
	// Sanitize input strings.
	$genome  = sanitize_POST("genome");
	$project = sanitize_POST("project");

	// Confirm if requested genome/project exists.
	if ($genome != "") {
		$genome_dir  = "users/".$user."/genomes/".$genome;
		if (!is_dir($genome_dir)) {
			// Genome doesn't exist, should never happen: Force logout.
			session_destroy();
			header('Location: .');
		}
		$target_dir = "../../users/".$user."/genomes/".$genome."/";
	} else if ($project != "") {
		$project_dir = "users/".$user."/projects/".$project;
		if (!is_dir($genome_dir)) {
			// Genome doesn't exist, should never happen: Force logout.
			session_destroy();
			header('Location: .');
		}
		$target_dir = "../../users/".$user."/projects/".$project."/";
	} else {
		// Genome and project strings sent at once, should never happen: Force logout.
		session_destroy();
		header('Location: .');
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
