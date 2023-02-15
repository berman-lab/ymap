<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	require_once 'POST_validation.php';
	ini_set('display_errors', 1);

	// Sanitize input strings.
	$user    = sanitize_POST("user");
	$pw_in   = stripHTML_POST("pw");

	// Validate login.
	validateLogin($user, $pw_in);

	// log-in success boolean.
	$login_success = false;

	// Load any delay from session; to prevent brute-force attacks by reloading page.
	if (isset($_SESSION['delay'])) {
		$delay = $_SESSION['delay'];
		if ($delay > 0) {
			sleep($delay);
		}
	}
//=========================================================
// Functions used to validate login credentials.
//---------------------------------------------------------
	function validateLogin($user, $pw_in){
		global $pepper;
		if (file_exists("users/".$user."/")) {
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

// Delay before page reload.
if ($login_success == false) {
	$_SESSION['delay'] = 5;
	echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
	echo "var intervalID = window.setInterval(reload_page, 5000);\n</script>\n";
} else {
	$_SESSION['delay'] = 0;
	echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
	echo "var intervalID = window.setInterval(reload_page, 500);\n</script>\n";
}
?>
