<?php
	session_start();
	require_once 'constants.php';
	ini_set('display_errors', 1);

	// Debugging information.
	print_r($GLOBALS);

	$tmp_file_name = $_FILES['Filedata']['tmp_name'];
	$ok = move_uploaded_file($tmp_file_name, '/path_to/new_filename');

	// This message will be passed to 'oncomplete' function
	echo $ok ? "OK" : "FAIL";

// Debugging output.
//    echo "<pre>    global vars    : "; print_r($GLOBALS); echo "</pre>";
?>
