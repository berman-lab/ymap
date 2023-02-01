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

	$user          = $_SESSION['user'];

	$user_dir = "users/".$user."/";
	if (is_dir($user_dir)) {
		// Requested user dir does exist: Delete user.
		rrmdir($user_dir);
		// Force logout.
		session_destroy();
		header('Location: .');
	} else {
		// User doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: .');
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
