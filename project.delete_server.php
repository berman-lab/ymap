<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	ini_set('display_errors', 1);

	$user    = $_POST["user"];
	$project = $_POST["project"];

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
