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
	$user   = $_SESSION['user'];

	// Sanitize input string.
	$hapmap = sanitize_POST("hapmap");

	// Confirm if requested hapmap exists.
	$dir = "users/".$user."/hapmaps/".$hapmap;
	if (is_dir($dir)) {
		// Requested hapmap dir does exist for logged in user: Delete installed hapmap.
		rrmdir($dir);
	} else {
		// Hapmap doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: .');
	}

	// Function for recursive rmdir, to clean out full hapmap directory.
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
