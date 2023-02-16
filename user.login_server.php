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
	$login_success = validateLogin($user, $pw_in);

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

			// Check if user account is locked.
			if (file_exists("users/".$user."/locked.txt")) {
				// Account is locked pending admin approval.
				echo "<font color=\"red\"><b>ERROR: Account is temporarily locked pending admin approval.</b></font><br>\n";
				echo "This may happen because account was newly registered or other issues.</br>\n";
				echo "(Main page will reload shortly...)<br>\n";
				echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
				echo "var intervalID = window.setInterval(reload_page, 5000);\n</script>\n";

				// Set login_success to 1 to prevent password failure delay.
				$login_success = 1;
			} else {
				// Account is active.

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
					$login_success = 1;
				} else {
					echo "<font color=\"red\"><b>ERROR: Input did not match a registered username & password combination.</b></font><br>\n";
					echo "(Main page will reload shortly...)<br>\n";
					echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
					echo "var intervalID = window.setInterval(reload_page, 5000);\n</script>\n";
				}
			}
		} else {
			//User doesn't exist
			echo "<font color=\"red\"><b>ERROR: Input did not match a registered username & password combination.</b></font><br>\n";
			echo "(Main page will reload shortly...)<br>\n";
			echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
			echo "var intervalID = window.setInterval(reload_page, 5000);\n</script>\n";
		}
		return $login_success;
	}

// Delay before page reload.
if ($login_success == 0) {
	$_SESSION['delay'] = 5;
	echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
	echo "var intervalID = window.setInterval(reload_page, 5000);\n</script>\n";
} else {
	$_SESSION['delay'] = 0;
	echo "<script type=\"text/javascript\">\nreload_page=function() {\n\tlocation.replace(\"panel.user.php\");\n}\n";
	echo "var intervalID = window.setInterval(reload_page, 500);\n</script>\n";
}
?>
