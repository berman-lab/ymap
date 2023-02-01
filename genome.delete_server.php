<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	require_once 'POST_validation.php';
	ini_set('display_errors', 1);

        // If the user is not logged on, redirect to login page.
        if(!isset($_SESSION['logged_on'])){
		session_destroy();
                header('Location: .');
        }

	// Load user string from session.
	$user    = $_SESSION['user'];

	// Sanitize input string.
	$genome  = sanitize_POST("newGenomeName");
	$dir     = "users/".$user."/genomes/".$genome;

	// Confirm if requested project exists.
	if (is_dir($project_dir)) {
		// Requested project dir does exist for logged in user: Delete installed project.
		rrmdir($dir);
		echo "COMPLETE";
	} else {
		// Project doesn't exist, should never happen. 
		echo "ERROR:".$user." doesn't own project.";
	}

	// Function for recursive rmdir, to clean out full genome directory.
	function rrmdir($dir) {
		if (is_dir($dir)) {
			$objects = scandir($dir);
			foreach ($objects as $object) {
				if ($object != "." && $object != "..") {
					if (filetype($dir."/".$object) == "dir") rrmdir($dir."/".$object); else unlink($dir."/".$object);
				}
			}
			reset($objects);
			rmdir($dir);
		}
	}
?>
