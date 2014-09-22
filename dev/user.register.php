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
		<title>Ymap | New User Registration</title>
	</head>
	<body>
		<b>New user registration.</b>
		<div class='tab' style='font-size:10pt'>
			The username doesn't need to correspond to your identity and passwords are encrypted, making them unavailable even to system administrators.<br>
			The true identity of the users can only be revealed if they choose to notify the administrators, which is generally not needed.<br><br>
		</div>
		<div style='font-size:10pt'>
		<form action="php/user.register_server.php" method="post">
			<label for="user">Username: </label>                                        <input type="text"     id="user"                      name="user"><br>
			<label for="pw1">Password: </label>                                         <input type="password" id="pwOrig"                    name="pwOrig">
			<label for="pw2">Retype Password: </label>                                  <input type="password" id="pwCopy"                    name="pwCopy"><br>
			<input type="submit" value="Register User">
			<input type="button" value="Cancel" onclick="location.replace('panel.user.php');">
		</form>
		</div>
	</body>
</html>
