<?php
	session_start();
	require_once 'constants.php';
	ini_set('display_errors', 1);

	$bad_chars = array("~","@","#","$","%","^","&","*","(",")","+","=","|","{","}","<",">","?",".",",","\\","/","'",'"',"[","]","!");
	$userOrig  = str_replace($bad_chars,"",trim(filter_input(INPUT_POST, "user",   FILTER_SANITIZE_STRING)));

	$primaryInvestigatorName  = filter_input(INPUT_POST, "primaryInvestigator_name",  FILTER_SANITIZE_STRING);
	$primaryInvestigatorEmail = filter_input(INPUT_POST, "primaryInvestigator_email", FILTER_SANITIZE_STRING);
	$researchInstitution      = filter_input(INPUT_POST, "researchInstitution",       FILTER_SANITIZE_STRING);
	$secondaryName            = filter_input(INPUT_POST, "secondary_name",            FILTER_SANITIZE_STRING);
	$secondaryEmail           = filter_input(INPUT_POST, "secondary_email",           FILTER_SANITIZE_STRING);
	$pwOrig                   = filter_input(INPUT_POST, "pwOrig",                    FILTER_SANITIZE_STRING);
	$pwCopy                   = filter_input(INPUT_POST, "pwCopy",                    FILTER_SANITIZE_STRING);

	$user = validateUser($userOrig);
	$pw   = validatePassword($pwOrig, $pwCopy);

	// User and Password both validated as correct
	if($user && $pw){
		createNewUser($user, $pw);
		createSecondaryInformationFile($user, $primaryInvestigatorName, $primaryInvestigatorEmail, $researchInstitution, $secondaryName, $secondaryEmail);
	}

//=========================================================
// Functions used to generate user account.
//---------------------------------------------------------
	function createNewUser($user, $pw){
		if (!doesUserDirectoryExist($user)) {
			$dir = "users/".$user;
			mkdir($dir);
			chmod($dir,0777);
			mkdir($dir."/projects/");    // initialize user projects dir.
			mkdir($dir."/genomes/");     // initialize user genomes dir.
			chmod($dir."/projects/",0777);
			chmod($dir."/genomes/", 0777);
			writePassword($user, $pw);
			$_SESSION['logged_on'] = 1;
			$_SESSION['user']      = $user;
			echo "<font color=\"green\"><b>SUCCESS: User account created.</b></font><br>";
			echo "(Main page will reload shortly...)\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 1000);\n</script>\n";
			return false;
		} else {
			// The user directory already exists!
			echo "<font color=\"red\"><b>ERROR: Invalid user name, try another.</b></font><br>";
			echo "(Main page will reload shortly...)\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"user.register.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 1000);\n</script>\n";
			return false;
		}
	}
	function writePassword($user, $pw){
		$pwFile = "users/".$user."/pw.txt";
		$fh     = fopen($pwFile, 'w');
		fwrite($fh, $pw);
		fclose($fh);
		chmod($pwFile, 0644);
	}
	function doesUserDirectoryExist($user){
		$dir = "users/".$user."/";
		return file_exists($dir);
	}
	function createSecondaryInformationFile($user, $primaryInvestigatorName, $primaryInvestigatorEmail, $researchInstitution, $secondaryName, $secondaryEmail){
		$secondaryInformationFile = "users/".$user."/info.txt";
		$fileHandle               = fopen($secondaryInformationFile, 'w');
		fwrite($fileHandle, "User: ".$user."\n");
		fwrite($fileHandle, "Primary Investigator Name: ".$primaryInvestigatorName."\n");
		fwrite($fileHandle, "Primary Investigator Email: ".$primaryInvestigatorEmail."\n");
		fwrite($fileHandle, "Research Institution: ".$researchInstitution."\n");
		fwrite($fileHandle, "Secondary Name: ".$secondaryName."\n");
		fwrite($fileHandle, "Secondary Email: ".$secondaryEmail."\n");
		fclose($fileHandle);
		chmod($secondaryInformationFile, 0644);
	}

//=========================================================
// Functions used to validate entered user name.
//---------------------------------------------------------
	function validateUser($user){
		$MIN_USER_LENGTH = 6;
		$MAX_USER_LENGTH = 24;
		// MIN LENGTH CHECK
		if (strlen($user) < $MIN_USER_LENGTH) {
			echo "<font color=\"red\"><b>ERROR: Your user name is too short, minimum is $MIN_USER_LENGTH.</b></font><br>";
			echo "(Main page will reload shortly...)\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"user.register.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 1000);\n";
			echo "</script>\n";
			return "";
		}
		// MAX LENGTH CHECK
		if (strlen($user) > $MAX_USER_LENGTH) {
			echo "<font color=\"red\"><b>ERROR: Your user name is too long, maximum is $MAX_USER_LENGTH.</b></font><br>";
			echo "(Main page will reload shortly...)\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"user.register.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 1000);\n";
			echo "</script>\n";
			return "";
		}
		//CHECK FOR NON ALPHANUMERIC CHARACTERS
		if (checkForAlphanumericCharacters($user)) {
			echo "<font color=\"red\"><b>ERROR: You have a non-alphanumeric character in your username.</b></font><br>";
			echo "(Main page will reload shortly...)\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"user.register.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 1000);\n";
			echo "</script>\n";
			return "";
		}
		//RETURN user
		return $user;
	}
	function checkForAlphanumericCharacters($string){
		return preg_match( "/^[a-zA-Z0-9]$/", $string);
	}

//=========================================================
// Function used to validate entered user password.
//---------------------------------------------------------
	function validatePassword($pwOrig, $pwCopy){
		$MIN_PASSWORD_LENGTH = 6;
		$MAX_PASSWORD_LENGTH = 24;
		// MIN LENGTH CHECK
		if (strlen($pwOrig) < $MIN_PASSWORD_LENGTH) {
			echo "<font color=\"red\"><b>ERROR: Your password is too short, minimum is $MIN_PASSWORD_LENGTH.</b></font><br>";
			echo "(Main page will reload shortly...)\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"user.register.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 1000);\n";
			echo "</script>\n";
			return "";
		}
		// MAX LENGTH CHECK
		if (strlen($pwOrig) > $MAX_PASSWORD_LENGTH) {
			echo "<font color=\"red\"><b>ERROR: Your password is too long, maximum is $MAX_PASSWORD_LENGTH.</b></font><br>";
			echo "(Main page will reload shortly...)\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"user.register.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 1000);\n";
			echo "</script>\n";
			return "";
		}
		// ORIG = COPY check
		if ($pwOrig != $pwCopy) {
			echo "<font color=\"red\"><b>ERROR: The passwords that you entered do not match.</b></font><br>";
			echo "(Main page will reload shortly...)\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"user.register.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 1000);\n";
			echo "</script>\n";
			return "";
		}
		return md5($pwOrig);
	}
?>
