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
Your account has been provided with administrator priviledges.<br><br>
You can now examine and respond to bug reports using the "Bug Reporting" tab.<br>
This "Admin" tab will have system troubleshooting tools and/or notes for administrators, but cannot be seen by normal user accounts.<br>
It may be helpful for you to make a second account to view the site as a normal user.<br>
</font>
<hr width="100%">

<font size='3'>
Eventually, this tab will provide account management options for the administrator(s).<br>
<ol>
	<li><b>View user</b>: See user project lists, for troubleshooting datasets and other purposes.</li>
	<li><b>Reset password</b>: Send the user an email with temporary password. Will require email addresses in user registration form (now done) and other system updates.</li>
	<li><b>User administration</b>: Change user/admin status of other users. Changing status currently requires server access.</li>
	<li><b>Delete user</b>: Delete inactive users. If a user has not been active for a still-to-be-determined length of time, the account should be deleted to clean up server space.
		In the future, the last date of activity will be displayed below.</li>
</ol>
</font>
<hr width="100%">

<font size='3'>
To assist in setting up Y<sub>MAP</sub> on another server, the default account contents can be downloaded as the following ZIP archives. The default account contents can be
regenerated using the pipeline tools and published datasets, but it takes a few days to do so.<br>
&nbsp; &nbsp; &nbsp; <a href="user_backup/default/genomes_backup.zip">genomes_backup.zip</a> [785 MB]<br>
&nbsp; &nbsp; &nbsp; <a href="user_backup/default/hapmaps_backup.zip">hapmaps_backup.zip</a> [1,817 KB]<br>
&nbsp; &nbsp; &nbsp; <a href="user_backup/default/projects_backup.zip">projects_backup.zip</a> [3.2 GB]<br>
The source files and directory structure of the pipeline website are managed through the GIT administrative tools I've mentioned. The files in these archives are too large for the GIT system
to manage effectively, so I had to put together this alternative.
</font>
<hr width="100%">

<font size='3'>
Future Updates to website:
<ol>
	<li>Add citation to "Help" tab when paper is published.</li>
	<li>Add BAM file download to "Manage Datasets" tab. A link at the right of existing download links, looking like "[BAM file]".</li>
</ol>
</font>
<hr width="100%">

The radio-buttons to the left of each user are intended to allow administrators to view that user's projects, similar to what is seen in the "Visualize Datasets" tab.
I have not yet gotten to implementing this function.
<table width="100%" cellpadding="0"><tr>
<td width="65%" valign="top">
	<?php
	//.---------------.
	//| User Accounts |
	//'---------------'
	if (isset($_SESSION['logged_on'])) {
		$userDir      = "users/";
		$userFolders  = array_diff(glob($userDir."*\/"), array('..', '.'));
		// Sort directories by date, newest first.
		array_multisort($userFolders, SORT_ASC, $userFolders);
		// Trim path from each folder string.
		foreach($userFolders as $key=>$folder) {   $userFolders[$key] = str_replace($userDir,"",$folder);   }
		$userCount = count($userFolders);

		echo "<table width='100%'>";
		echo "<tr><td width='20%'><font size='2'><b>User Account</b></font></td><td><font size='2'><b>Last Active</b></font></td></tr>";
		foreach($userFolders as $key=>$user) {
			echo "<tr><td><span id='project_label_".$key."' style='color:#000000;'>\n\t\t";
			echo "<font size='2'>".($key+1).".";
			echo "<input  id='show_".$key."' type='checkbox'\">";
			echo "\n\t\t".$user."</font></span>\n\t\t";
			echo "<span id='u2_".$user."_delete'></span><span id='u_".$user."_type'></span>\n\t\t";
			echo "<br>\n\t\t";
			echo "<div id='frameContainer.u1_".$key."'></div>";
			echo "</td><td>---</td></tr>";
		}
		echo "</table>";
	} else {
		$userCount = 0;
	}
	?>
</td><td width="35%" valign="top">
</td></tr></table>
</html>
