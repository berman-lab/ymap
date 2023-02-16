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
This "Admin" tab will have system troubleshooting tools and/or notes for administrators, but cannot be seen by normal user accounts.
It may be helpful for you to make a second account to view the site as a normal user.<br>
</font>
<hr width="100%">

Newly registered locked accounts, pending approval.
<table width="100%" cellpadding="0"><tr>
<td width="65%" valign="top">
<script type="text/javascript" src="js/jquery-3.6.3.js"></script>
<script type="text/javascript" src="js/jquery.form.js"></script>
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
		echo "<tr><td width='20%'><font size='2'><b>User Account</b></font></td><td><font size='2'><b>Approval Needed</b></font></td></tr>";
		foreach($userFolders as $key=>$user) {
			echo "<tr><td><span id='project_label_".$key."' style='color:#000000;'>\n\t\t";
			echo "<font size='2'>".($key+1).".";

			echo "\n\t\t".$user."</font></span>\n\t\t";
			echo "<span id='u2_".$user."_delete'></span><span id='u_".$user."_type'></span>\n\t\t";
			echo "<br>\n\t\t";
			echo "<div id='frameContainer.u1_".$key."'></div>";
			echo "</td><td>";

			if (file_exists("users/".$user."/locked.txt")) {
				echo "\n<input type='button' value='Approve' onclick=\"key = '$key'; console.log(key); $.ajax({url:'admin.approve_server.php',type:'post',data:{key:key},success:function(answer){console.log(answer);}});location.replace('panel.admin.php');\">\n";
			}

			echo "</td></tr>";
		}
		echo "</table>";
	} else {
		$userCount = 0;
	}
	?>
</td><td width="35%" valign="top">
</td></tr></table>
</html>
