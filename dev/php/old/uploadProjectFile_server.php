<?php
	session_start();
	require_once 'constants.php';
	ini_set('display_errors', 1);

	//echo "Success!\n";
	if ($_FILES["fileToBeUploaded"]["error"] > 0){
		echo "ERROR";
	} else {
		$temporaryFileName = $_FILES["fileToBeUploaded"]["tmp_name"];
		$project           = str_replace("\\",".",str_replace("/",".",str_replace(" ","_",$_POST["upload_projectName"])));
		$filenameOnServer  = $_POST["upload_filenameOnServer"];
		$user              = $_SESSION['user'];
		$successReturn     = array($project, $filenameOnServer);

		if(move_uploaded_file($temporaryFileName, $GLOBALS['directory']."users/".$user."/projects/".$project."/".$filenameOnServer)){
			//Move file was successful
			chmod($GLOBALS['directory']."users/".$user."/projects/".$project."/".$filenameOnServer, 0777);
			echo json_encode($successReturn);
		} else {
			//Move file was unsuccessful
			echo "ERROR_FILE_MOVE";
		}
	}
?>
