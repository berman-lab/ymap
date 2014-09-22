<?php
	session_start();
	require_once 'constants.php';
	ini_set('display_errors', 1);

    print_r($GLOBALS);

	$temporaryFileName = $_FILES["fileToBeUploaded"]["tmp_name"];
	$genome            =  str_replace("\\",".",str_replace("/",".",str_replace(" ","_",$_POST["upload_genomeName"])));
	$filenameOnServer  = $_FILES["fileToBeUploaded"]['name']; // $_POST["upload_filenameOnServer"];
	$user              = $_SESSION['user'];
	$successReturn     = array($genome, $filenameOnServer);

	echo "<pre>Debugging Output:</pre>";
	echo "<pre>    temp file name      : ".$temporaryFileName."</pre>";
	echo "<pre>    genome              : ".$genome."</pre>";
	echo "<pre>    file name on server : "; print_r($filenameOnServer); echo "</pre>";
	echo "<pre>    user                : ".$user."</pre>";
	echo "<pre>    sucess return       : "; print_r($successReturn); echo"</pre>";

	if ($_FILES["fileToBeUploaded"]["error"] > 0){
		echo "ERROR";
	} else {
		// call "move_uploaded_file()" with details of file for transport.
		if(move_uploaded_file($temporaryFileName, $GLOBALS['directory']."users/".$user."/genomes/".$genome."/".$filenameOnServer)){
			//Move file was successful
			chmod($GLOBALS['directory']."users/".$user."/genomes/".$genome."/".$filenameOnServer, 0777);
			echo json_encode($successReturn);
		} else {
			//Move file was unsuccessful
			echo "ERROR_FILE_MOVE : File transfer was not completed.";
		}
	}
?>
