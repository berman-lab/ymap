<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';

	ini_set('display_errors', 1);

	$bad_chars = array("~","@","#","$","%","^","&","*","(",")","+","=","|","{","}","<",">","?",".",",","\\","/","'",'"',"[","]","!");
	$user      = str_replace($bad_chars,"",trim(filter_input(INPUT_POST, "user",   FILTER_SANITIZE_STRING)));

	if(isset($_SESSION['logged_on']) && $user == $_SESSION['user']){
		// User confirmed, can delete user.
		$dir = "users/".$user."/";
		rrmdir($dir);
		session_unset();
		echo "COMPLETE";
	} else {
		echo "ERROR";
	}

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
