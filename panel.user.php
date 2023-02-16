<?php
	session_start();
	if (isset($_SESSION['logged_on'])) { $user = $_SESSION['user']; } else { $user = 'default'; }
?>
<style type="text/css">
html * {
	font-family: arial !important;
}
</style>
<span id="firefox_error_span"></span>
<script type="text/javascript">
if (navigator.userAgent.toLowerCase().indexOf('firefox') > -1) {
	// Do Firefox-related activities
	var FFerror       = document.getElementById("firefox_error_span");
	ErrorString       = '<font color="red"><b>';
	ErrorString      += 'Some features of this website require a web-browser based on the Blink (Chrome, Opera, etc.) or WebKit (Safari, etc.) rendering engines.<br><br>';
	ErrorString      += 'Firefox is based on the Gecko rendering engine. YMAP has an error when rendered with this engine, resulting in datasets not processing after data upload.';
	ErrorString      += '</b></font><br><br>';
	FFerror.innerHTML = ErrorString;
}
</script>

<font size='3'>Log into a preexisting user account or create a new user account.</font><br>
<?php
if ((isset($_SESSION['delay'])) && !(isset($_SESSION['logged_on']))) {
	$delay = $_SESSION['delay'];
	if ($delay != 0) {
		echo "<font size='2' color='Red'>(There will be a short delay afer hitting 'Log In' button due to prior log in failure.)</font><br><br>";
	}
}
echo "<br>";

if (isset($_SESSION['logged_on'])) {
	echo "User '<b>".$user."</b>' logged in. \n";
	// provide logout button.
	echo "<button type='button' onclick=\"window.location.href='user.logout_server.php'\">Logout</button>\n";
	// provide delete-user button.
	echo "<span id='u_".$user."_delete'></span>\n";
	echo "<button type='button' onclick=\"window.location.href='user.delete.php'\">Delete User.</button>\n";
	echo "<br><br>\n";
	echo "<font size='2'>\n\t";
	echo "You can navigate through the above menu and show/close projects while new datafiles are uploading.<br>\n\t";
	echo "A page reload or project/genome/hapmap creation/deletion, however, will interrupt file transfer.<br><br>\n\t";
	echo "Depending on system load, tasks may take an hour or more to complete after data upload is complete.<br><br>\n\t";
	echo "Reload page and select 'projects' tab to check for newly completed projects.\n";
	echo "</font>\n";

	if (isset($_SESSION['reload_once'])) {
		echo "<script type=\"text/javascript\"> parent.location.reload(); </script>\n";
		unset($_SESSION['reload_once']);
	}
} else {
	echo "<script type=\"text/javascript\">\n\t\n\t</script>\n\t";
	echo "<form action='user.login_server.php' method='post'>\n\t";
	echo "<label for='user'>Username: </label><input type='text' id='user' name='user'><br>\n\t";
	echo "<label for='pw'>Password: </label><input type='password' id='pw' name='pw'><br>\n\t";
	echo "<button type='submit'>Log In</button>\n\t";
	echo "</form>\n\t";
	echo "<font size='2'>";
	echo "If you don't have a user account, you may make one by clicking below.<br>Your account will initially remain locked until it has been approved by an admin.<br>\n\t";
	echo "<button type='button' onclick=\"window.location.href='user.register.php'\">Register new user.</button>\n\t";
	echo "</font>\n\t";
	$_SESSION['reload_once'] = "true";
}
?>
