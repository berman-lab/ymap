<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	ini_set('display_errors', 1);

	$bad_chars = array("~","@","#","$","%","^","&","*","(",")","+","=","|","{","}","<",">","?",".",",","\\","/","'",'"',"[","]","!");
	$user      = str_replace($bad_chars,"",trim(filter_input(INPUT_POST, "user",   FILTER_SANITIZE_STRING)));

	$pw_in   = filter_input(INPUT_POST, "pw", FILTER_SANITIZE_STRING);

	$currentPath = getcwd();
	$user        = validateUser($user_in);
	$pw          = validatePassword($pw_in);

	validateLogin($user, $pw);

//=========================================================
// Functions used to validate login credentials.
//---------------------------------------------------------
	function validateLogin($user, $pw){
		global $currentPath;
		if(doesUserDirectoryExist($user)){
			// User Exists
			$userDir = "users/".$user."/";
			$pwFile  = "users/".$user."/pw.txt";

			if (file_exists($userDir)) {
				if(file_get_contents($pwFile) === $pw) {
					$_SESSION['logged_on'] = 1;
					$_SESSION['user']      = $user;
					echo "<font color=\"green\"><b>SUCCESS: User is now logged in.</b></font><br>\n";
					echo "(Main page will reload shortly...)<br>\n";
					echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
					echo "var intervalID = window.setInterval(reload_page, 1000);\n</script>\n";
				} else {
					echo "<font color=\"red\"><b>ERROR: Input did not match a registed username & password combination.</b></font><br>\n";
					echo "(Main page will reload shortly...)<br>\n";
					echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
					echo "var intervalID = window.setInterval(reload_page, 1000);\n</script>\n";
				}
			} else {
				echo "<font color=\"red\"><b>ERROR: Input did not match a registed username & password combination.</b></font><br>\n";
				echo "(Main page will reload shortly...)<br>\n";
				echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
				echo "var intervalID = window.setInterval(reload_page, 1000);\n</script>\n";
			}
		} else {
			//User doesn't exist
			echo "<font color=\"red\"><b>ERROR: Input did not match a registered username & password combination.</b></font><br>\n";
			echo "(Main page will reload shortly...)<br>\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 1000);\n</script>\n";
		}
	}
	function doesUserDirectoryExist($user){
		$dir = "users/".$user."/";
		return file_exists($dir);
	}

//=========================================================
// Functions used to validate entered user name.
//---------------------------------------------------------
	function validateUser($user){
		$MIN_USER_LENGTH = 6;
		$MAX_USER_LENGTH = 24;
		// MIN LENGTH CHECK
		if(strlen($user) < $MIN_USER_LENGTH){
			echo "<font color=\"red\"><b>ERROR: Usernames must be at least $MIN_USER_LENGTH characters long.</b></font><br>\n";
			echo "(Main page will reload shortly...)<br>\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 1000);\n</script>\n";
			return "";
		}
		// MAX LENGTH CHECK
		if(strlen($user) > $MAX_USER_LENGTH){
			echo "<font color=\"red\"><b>ERROR: Usernames must be at most $MAX_USER_LENGTH characters long.</b></font><br>\n";
			echo "(Main page will reload shortly...)<br>\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 1000);\n</script>\n";
			return "";
		}
		//CHECK FOR NON ALPHANUMERIC CHARACTERS
		if(checkForAlphanumericCharacters($user)){
			echo "<font color=\"red\"><b>ERROR: Your username contains non-alphanumeric characters.</b></font><br>\n";
			echo "(Main page will reload shortly...)<br>\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 1000);\n</script>\n";
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
	function validatePassword($pw){
		$MIN_PASSWORD_LENGTH = 6;
		$MAX_PASSWORD_LENGTH = 24;
		// MIN LENGTH CHECK
		if(strlen($pw) < $MIN_PASSWORD_LENGTH){
			echo "<font color=\"red\"><b>ERROR: Passwords must be at least $MIN_PASSWORD_LENGTH.</b></font><br>\n";
			echo "(Main page will reload shortly...)<br>\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 1000);\n</script>\n";
			return "";
		}
		// MAX LENGTH CHECK
		if(strlen($pw) > $MAX_PASSWORD_LENGTH){
			echo "<font color=\"red\"><b>ERROR: Passwords must be at most $MAX_PASSWORD_LENGTH.</b></font><br>\n";
			echo "(Main page will reload shortly...)<br>\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 1000);\n</script>\n";
			return "";
		}
		return md5($pw);
	}
?>
<script type="text/javascript">
	parent.location.reload();
</script>
