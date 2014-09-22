<?php
	session_start();
	require_once 'php/constants.php';
?>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
        "http://www.w3.org/TR/html4/loose.dtd">

<html lang="en">
	<head>
		<style type="text/css">
			body {font-family:  arial;}
			.tab { margin-left: 1cm;  }
		</style>
		<meta http-equiv="content-type" content="text/html; charset=utf-8">
		<title>Ymap | Delete User</title>
		<!-- Needed for deletion functions.--!>
		<script type="text/javascript" src="js/jquery-1.7.2.js"></script>
		<script type="text/javascript" src="js/jquery.form.js"></script>
		<script type="text/javascript">
		user = '<?php echo $_SESSION['user']; ?>';
		</script>
	</head>
	<body>
		<b>Delete User.</b>
		<div class='tab' style='font-size:10pt'>
			This action will <b>PERMANENTLY</b> delete the '<i><?php echo $_SESSION['user']; ?></i>' user account and all associated files.<br>
			There are no backups to restore from.<br>
			Are you absolutely sure?<br>
			<div style='font-size:10pt'>
			<form action="$.ajax({url : 'php/user.delete_server.php', type : 'post', data : {user : user}, success : function(answer){if(answer == 'COMPLETE'){window.location.href=window.location.href;}} });" method="post">
				<input type="button" value="Yes! Delete." onclick="console.log('@user.delete.php@ user = '+user); $.ajax({url : 'php/user.delete_server.php', type : 'post', data : {user : user}, success : function(answer){if(answer == 'COMPLETE'){window.location.href=window.location.href;}} }); location.replace('panel.user.php');">
				<input type="button" value="No! Cancel."  onclick="location.replace('panel.user.php');">
			</form>
			</div>
		</div>
	</body>
</html>
