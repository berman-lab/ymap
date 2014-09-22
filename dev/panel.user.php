<?php
	session_start();
	if (isset($_SESSION['logged_on'])) { $user = $_SESSION['user']; } else { $user = 'default'; }
?>
<style type="text/css">
html * {
	font-family: arial !important;
}
</style>
<font size='3'>Log into a preexisting user account or create a new user account.</font><br><br>
<?php
if (isset($_SESSION['logged_on'])) {
	echo "User '<b>".$user."</b>' logged in. \n";
	// provide logout button.
	echo "<button type='button' onclick=\"window.location.href='php/user.logout_server.php'\">Logout</button>\n";
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
} else {
	echo "<script type=\"text/javascript\">\n\t\n\t</script>\n\t";
	echo "<form action='php/user.login_server.php' method='post'>\n\t";
	echo "<label for='user'>Username: </label><input type='text' id='user' name='user'><br>\n\t";
	echo "<label for='pw'>Password: </label><input type='password' id='pw' name='pw'><br>\n\t";
	echo "<button type='submit' onclick=\"parent.update_projectsShown_after_logout();\">Login</button>\n\t";
	echo "</form>\n\t";
	echo "<font size='2'>";
	echo "If you don't have a user account, you may make one by clicking below.<br>\n\t";
	echo "<button type='button' onclick=\"window.location.href='user.register.php'\">Register new user.</button>\n\t";
	echo "</font>\n\t";
}
?>
