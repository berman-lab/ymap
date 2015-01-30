<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.
	// Exec random.pl which will return true or false based on whether there is a problem
	$exec_str = "scripts/random.pl 2>&1";
	exec($exec_str, $exec_output, $exec_return);

	if(!$exec_return){
		echo json_encode($exec_output);
	} else {
		echo "FAILURE";	
	}


?>
