<?php
	session_start();
	if (isset($_SESSION['logged_on'])) { $user = $_SESSION['user']; } else { $user = 'default'; }
?>
<html style="background: #FFDDDD;">
<style type="text/css">
	html * {
		font-family: arial !important;
	}
</style>
<font size='3'>
Eventually, this tab will provide the user account management options for the administrator.<br>
<ol>
	<li><b>View user</b>: See site from perspective of user's account, for troubleshooting datasets and other purposes.</li>
	<li><b>Reset password</b>: Send the an email with temporary password. Will require email addresses in user registration.</li>
</ol>
Currently, it simply shows a list of registered user accounts.</font><br><br>
<table width="100%" cellpadding="0"><tr>
<td width="65%" valign="top">
	<?php
	//.---------------.
	//| User Accounts |
	//'---------------'
	if (isset($_SESSION['logged_on'])) {
		$userDir      = "users/";
		$userFolders   = array_diff(glob($userDir."*"), array('..', '.'));
		// Sort directories by date, newest first.
		array_multisort($userFolders, SORT_ASC, $userFolders);
		// Trim path from each folder string.
		foreach($userFolders as $key=>$folder) {   $userFolders[$key] = str_replace($userDir,"",$folder);   }
		$userCount = count($userFolders);

		echo "<b><font size='2'>Registered User accounts:</font></b>\n\t\t";
		echo "<br>\n\t\t";
		foreach($userFolders as $key=>$user) {
			echo "<div class='tab'>";
			echo "<span id='project_label_".$key."' style='color:#000000;'>\n\t\t";
			echo "<font size='2'>".($key+1).".";
			echo "<input  id='show_".$key."' type='checkbox'\">";
			echo "\n\t\t".$user."</font></span>\n\t\t";
			echo "<span id='u2_".$user."_delete'></span><span id='u_".$user."_type'></span>\n\t\t";
			echo "<br>\n\t\t";
			echo "<div id='frameContainer.u1_".$key."'></div>";
			echo "</div>";
		}
	} else {
		$userCount = 0;
	}
	?>
</td><td width="35%" valign="top">
</td></tr></table>
</html>
