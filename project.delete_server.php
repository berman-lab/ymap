<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	ini_set('display_errors', 1);

        // If the user is not logged on, redirect to login page.
        if(!isset($_SESSION['logged_on'])){
                header('Location: user.login.php');
        }

	$user    = $_SESSION['user'];
        $project = trim(filter_input(INPUT_POST, "project", FILTER_SANITIZE_STRING));  // removing unwanted characters

	if($user == $_SESSION['user']){
		// User confirmed, can delete project
		$dir = "users/".$user."/projects/".$project;
		rrmdir($dir);
		echo "COMPLETE";
	} else {
		echo "ERROR:".$user;
	}

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
