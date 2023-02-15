<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	require_once 'POST_validation.php';
	ini_set('display_errors', 1);

	// Sanitize input strings.
	$user_in = sanitize_POST("user");
	$pw_in   = stripHTML_POST("pw");

	// Validate user and password inputs.
	$user    = validateUser($user_in);
	$pw_hash = validatePassword($pw_in);

	// Validate login.
	validateLogin($user, $pw_in);

	// log-in success boolean.
	$login_success = false;
//=========================================================
// Functions used to validate login credentials.
//---------------------------------------------------------
	function validateLogin($user, $pw_in){
		global $pepper;
		if (doesUserDirectoryExist($user)) {
			// User exists, so we check password.

			// Load stored password hash.
			$pwFile         = "users/".$user."/pw.txt";
			$pw_stored_hash = file_get_contents($pwFile);

			// Compare peppered input password to stored hash.
			$checked = password_verify($pw_in.$pepper, $pw_stored_hash);
			if ($checked) {
				$_SESSION['logged_on'] = 1;
				$_SESSION['user']      = $user;
				echo "<font color=\"green\"><b>SUCCESS: User is now logged in.</b></font><br>\n";
				echo "(Main page will reload shortly...)<br>\n";
				echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
				echo "var intervalID = window.setInterval(reload_page, 1000);\n</script>\n";
				$login_success = true;
			} else {
				echo "<font color=\"red\"><b>ERROR: Input did not match a registered username & password combination.</b></font><br>\n";
				echo "(Main page will reload shortly...)<br>\n";
				echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
				echo "var intervalID = window.setInterval(reload_page, 5000);\n</script>\n";
			}
		} else {
			//User doesn't exist
			echo "<font color=\"red\"><b>ERROR: Input did not match a registered username & password combination.</b></font><br>\n";
			echo "(Main page will reload shortly...)<br>\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 5000);\n</script>\n";
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
			echo "var intervalID = window.setInterval(reload_page, 5000);\n</script>\n";
			return "";
		}
		// MAX LENGTH CHECK
		if(strlen($user) > $MAX_USER_LENGTH){
			echo "<font color=\"red\"><b>ERROR: Usernames must be at most $MAX_USER_LENGTH characters long.</b></font><br>\n";
			echo "(Main page will reload shortly...)<br>\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 5000);\n</script>\n";
			return "";
		}
		//CHECK FOR NON ALPHANUMERIC CHARACTERS
		if(checkForAlphanumericCharacters($user)){
			echo "<font color=\"red\"><b>ERROR: Your username contains non-alphanumeric characters.</b></font><br>\n";
			echo "(Main page will reload shortly...)<br>\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 5000);\n</script>\n";
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
		global $pepper;
		$MIN_PASSWORD_LENGTH = 6;
		$MAX_PASSWORD_LENGTH = 24;
		// MIN LENGTH CHECK
		if(strlen($pw) < $MIN_PASSWORD_LENGTH){
			echo "<font color=\"red\"><b>ERROR: Passwords must be at least $MIN_PASSWORD_LENGTH.</b></font><br>\n";
			echo "(Main page will reload shortly...)<br>\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 5000);\n</script>\n";
			return "";
		}
		// MAX LENGTH CHECK
		if(strlen($pw) > $MAX_PASSWORD_LENGTH){
			echo "<font color=\"red\"><b>ERROR: Passwords must be at most $MAX_PASSWORD_LENGTH.</b></font><br>\n";
			echo "(Main page will reload shortly...)<br>\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 5000);\n</script>\n";
			return "";
		}

		// return modern, random-salted-peppered-hash.
		$peppered_pw = $pw.$pepper;
		return password_hash($peppered_pw, PASSWORD_DEFAULT, ['cost' => 10]);
	}

// Delay before page reload.
if ($login_success == false) {
	echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
	echo "var intervalID = window.setInterval(reload_page, 5000);\n</script>\n";
} else {
	echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
	echo "var intervalID = window.setInterval(reload_page, 500);\n</script>\n";
}
?>
