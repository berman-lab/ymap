<?php
	session_start();
	require_once 'constants.php';
	require_once 'POST_validation.php';
	require_once 'SecureNewDirectory.php';
	ini_set('display_errors', 1);

	// validate POST input.
	$primaryInvestigatorName  = sanitizeName_POST("primaryInvestigator_name");
	$primaryInvestigatorEmail = sanitizeEmail_POST("primaryInvestigator_email");

	// Not used, but for later use maybe.
	$researchInstitution      = sanitizeName_POST("researchInstitution");
	$secondaryName            = sanitizeName_POST("secondary_name");
	$secondaryEmail           = sanitizeEmail_POST("secondary_email");

	// user and password strings are validated in remainint code block.
	$userOrig   = filter_input(INPUT_POST, "user",        FILTER_SANITIZE_STRING);
	$pwOrig     = filter_input(INPUT_POST, "pwOrig",      FILTER_SANITIZE_STRING);
	$pwCopy     = filter_input(INPUT_POST, "pwCopy",      FILTER_SANITIZE_STRING);

	$user       = validateUser($userOrig);
	if ($user) { // User validated.
		$pw = validatePassword($pwOrig, $pwCopy);
	}
	if ($user && $pw) { // User and Password both validated.
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
			secureNewDirectory($dir);
			chmod($dir,0777);
			mkdir($dir."/projects/");    // initialize user projects dir.
			secureNewDirectory($dir."/projects/");
			mkdir($dir."/genomes/");     // initialize user genomes dir.
			secureNewDirectory($dir."/genomes/");
			mkdir($dir."/hapmaps/");     // initialize user hapmaps dir.
			secureNewDirectory($dir."/hapmaps/");
			chmod($dir."/projects/",0777);
			chmod($dir."/genomes/", 0777);
			writePassword($user, $pw);
			$_SESSION['logged_on'] = 1;
			$_SESSION['user']      = $user;
			echo "<font color=\"green\"><b>SUCCESS: User account created.</b></font><br>";
			echo "(Main page will reload shortly...)\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 500);\n</script>\n";
			return false;
		} else {
			// The user directory already exists!
			echo "<font color=\"red\"><b>ERROR: Invalid user name, try another.</b></font><br>";
			echo "(Main page will reload shortly...)\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"user.register.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 5000);\n</script>\n";
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
			echo "var intervalID = window.setInterval(reload_page, 5000);\n";
			echo "</script>\n";
			return "";
		}
		// MAX LENGTH CHECK
		if (strlen($user) > $MAX_USER_LENGTH) {
			echo "<font color=\"red\"><b>ERROR: Your user name is too long, maximum is $MAX_USER_LENGTH.</b></font><br>";
			echo "(Main page will reload shortly...)\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"user.register.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 5000);\n";
			echo "</script>\n";
			return "";
		}
		//CHECK FOR NON ALPHANUMERIC CHARACTERS
		if (preg_match('/[^a-zA-Z0-9]+/', $user)) {
			echo "<font color=\"red\"><b>ERROR: You have a non-alphanumeric character in your username.</b></font><br>";
			echo "(Main page will reload shortly...)\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"user.register.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 5000);\n";
			echo "</script>\n";
			return "";
		}
		//RETURN user
		return $user;
	}

//=========================================================
// Function used to validate entered user password.
//---------------------------------------------------------
	function validatePassword($pwOrig, $pwCopy){
		global $pepper, $user, $primaryInvestigatorName, $primaryInvestigatorEmail;
		//----------------------------------------------------------------
		// Checks to see if password was entered the same for both tries.
		//................................................................
		if ($pwOrig != $pwCopy) {
			echo "<font color=\"red\"><b>ERROR: The passwords that you entered do not match.</b></font><br>";
			echo "(Main page will reload shortly...)\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"user.register.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 5000);\n";
			echo "</script>\n";
			return "";
		}

		//-----------------------------------------------
		// Checks password validity by various criteria.
		//...............................................

		// MIN LENGTH CHECK
		$MIN_PASSWORD_LENGTH = 8;
		if (strlen($pwOrig) < $MIN_PASSWORD_LENGTH) {
			echo "<font color=\"red\"><b>ERROR: Your password is too short, minimum is $MIN_PASSWORD_LENGTH.</b></font><br>";
			echo "(Main page will reload shortly...)\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"user.register.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 5000);\n";
			echo "</script>\n";
			return "";
		}

		// MAX LENGTH CHECK
		$MAX_PASSWORD_LENGTH = 64;
		if (strlen($pwOrig) > $MAX_PASSWORD_LENGTH) {
			echo "<font color=\"red\"><b>ERROR: Your password is too long, maximum is $MAX_PASSWORD_LENGTH.</b></font><br>";
			echo "(Main page will reload shortly...)\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"user.register.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 5000);\n";
			echo "</script>\n";
			return "";
		}

		// Check to see if username is included in password.
		if (stristr(strtolower($pwOrig),strtolower($user))) {
			echo "<font color=\"red\"><b>ERROR: Your password cannot include your user name.</b></font><br>";
			echo "(Main page will reload shortly...)\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"user.register.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 5000);\n";
			echo "</script>\n";
			return "";
		}

		// Check to see if username tokens are included in password. Tokens split by any of: ",.-_ #\t"
		$name_tokens = preg_split( "/[,.-_ #;\t]/", $primaryInvestigatorName);
		for ($i=0; $i < sizeof($name_tokens); $i++) {
			if (strlen($name_tokens[$i]) > 3) {
				if (stristr(strtolower($pwOrig),strtolower($name_tokens[$i]))) {
					echo "<font color=\"red\"><b>ERROR: Your password cannot include parts of your name.</b></font><br>";
					echo "(Main page will reload shortly...)\n";
					echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"user.register.php\");\n}\n";
					echo "var intervalID = window.setInterval(reload_page, 5000);\n";
					echo "</script>\n";
					return "";
				}
			}
		}

		// Check to see if email address tokens are included in password. Tokens split by any of: ",.-_ #;\t@"
		$email_tokens = preg_split( "/[,.-_ #;\t@]/", $primaryInvestigatorEmail);
		for ($i=0; $i < sizeof($email_tokens); $i++) {
			if (strlen($email_tokens[$i]) > 3) {
				if (stristr(strtolower($pwOrig),strtolower($email_tokens[$i]))) {
					echo "<font color=\"red\"><b>ERROR: Your password cannot include parts of your email.</b></font><br>";
					echo "(Main page will reload shortly...)\n";
					echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"user.register.php\");\n}\n";
					echo "var intervalID = window.setInterval(reload_page, 5000);\n";
					echo "</script>\n";
					return "";
				}
			}
		}

		// Password must include characters in at least three of the following categories.
		// 1. Uppercase letters of European languages.
		// 2. Lowercase letters of European languages.
		// 3. Base 10 digis.
		// 4. Non-alphanumeri characters: ~!@#$%^&*_-+=`|\(){}[]:;"'<>,.?/
		//	Currency symbols such as the Euro or British Pound aren't counted as special characters for this policy setting.
		$char_req_1 = 0;
		$char_req_2 = 0;
		$char_req_3 = 0;
		$char_req_4 = 0;
		if (preg_match('/[A-Z]/', $pwOrig)) {
			// There is at least one uppercase letter.
			$char_req_1 = 1;
		}
		if (preg_match('/[a-z]/', $pwOrig)) {
			// There is at least one lowercase letter.
			$char_req_2 = 1;
		}
		if (preg_match('/[0-9]/', $pwOrig)) {
			// There is at least one numeral.
			$char_req_3 = 1;
		}
		if (preg_match('/[^a-zA-Z0-9]+/', $pwOrig)) {
			// There is at least one special character.
			$char_req_4 = 1;
		}
		if (($char_req_1+$char_req_2+$char_req_3+$char_req_4) < 3) {
			echo "<font color=\"red\"><b>ERROR: Your password must include characters from at least three of the following categories.</b></font><br>";
			echo "1. Capital letters.\n";
			echo "2. Lower case letters.\n";
			echo "3. Numerals.\n";
			echo "4. Special characters.\n\n";
			echo "(Main page will reload shortly...)\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"user.register.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 5000);\n";
			echo "</script>\n";
			return "";
		}

		// If the password has survived all the above checks, generate a random-salted-peppered-hash.
		$_SESSION['delay'] = 0;
		$peppered_pw = $pwOrig.$pepper;
		return password_hash($peppered_pw, PASSWORD_DEFAULT, ['cost' => 10]);
	}
?>
