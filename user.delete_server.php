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
	$dir     = "users/".$user."/";

	// Confirm if requested project exists.
	if (is_dir($project_dir)) {
		// Requested user does exist: Delete installed user.
		rrmdir($dir);
		echo "COMPLETE";
	} else {
		// Project doesn't exist, should never happen.
		echo "ERROR:".$user." doesn't exist.";
	}

	// recursive rmdir function.
	function rrmdir($dir) {
		if (is_dir($dir)) {
			$objects = scandir($dir);
			foreach ($objects as $object) {
				if ($object != "." && $object != "..") {
					if (filetype($dir."/".$object) == "dir") {
						rrmdir($dir."/".$object);
					} else {
						unlink($dir."/".$object);
					}
				}
			}
			reset($objects);
			rmdir($dir);
		}
	}
?>
