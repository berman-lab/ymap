<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	ini_set('display_errors', 1);

        // If the user is not logged on, redirect to login page.
        if(!isset($_SESSION['logged_on'])){
                header('Location: user.login.php');
        }

	$user = $_SESSION['user'];
//	$submitter   = filter_input(INPUT_POST, "submitter",   FILTER_SANITIZE_STRING);
	$description = filter_input(INPUT_POST, "description", FILTER_SANITIZE_STRING);

	$admin_email = '';

	if($user != "" && $description != ""){
		// Order of replacement
		$order   = array("\r\n", "\n", "\r");
		$replace = '<br>';

		// Processes \r\n's first so they aren't converted twice.
		$description = str_replace($order, $replace, $description);
		$description = str_replace("|", ":", $description);

		$bugtrackFile = "bugtracker/bugs.txt";

		// Write the new bug into the tracking file.
		$bugFileHandle = fopen($bugtrackFile, 'a');
		fwrite($bugFileHandle, $user."|".$description."| | |0\n");
		fclose($bugFileHandle);
		chmod("bugtracker/bugs.txt", 0777);

		// Email system administrator about new bug/reature comment.
		$email_subject     = "Ymap | Bug Report or Feature Request";
		$email_message     = "User account '".$user."' has added the following comment to the bug-tracking system.\n\n'".$description."'.\n-Ymap";
		$email_from        = "From: Ymap";
		mail($admin_email, $email_subject, $email_message);
	}
?>
