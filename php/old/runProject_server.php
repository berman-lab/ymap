<?php
	session_start();
	require_once 'constants.php';
	ini_set('display_errors', 1);

	$project = filter_input(INPUT_POST, "project", FILTER_SANITIZE_STRING);
	$ploidy  = filter_input(INPUT_POST, "ploidy",  FILTER_SANITIZE_NUMBER_FLOAT);
	$user    = $_SESSION['user'];

	//Validate ploidy it should be a number greater than 0
	if(!is_numberic($ploidy) || $ploidy < 0){
		$ploidy = 2;
	}
	if($project != "" && $user != ""){
		$input_directory = "/heap/hapmap/users/".$user."/projects/".$project;
		// Make sure everything exists properly then run if it does.
		if(file_exists($input_directory) && is_dir($input_directory)){
			exec('ssh-agent '.$GLOBALS['directory'].'scripts/runProject.sh '.escapeshellarg($input_directory).' '.escapeshellarg($project), $exec_output0, $exec_return0);
			echo "EXEC0: ".$exec_return0."\n";
			var_dump($exec_output0);
		} else {
			echo "VERY BAD!";
		}
	} else {
		echo "BAD!";
	}
?>
