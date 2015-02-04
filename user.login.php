<?php
	session_start();
//test
?>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
        "http://www.w3.org/TR/html4/loose.dtd">
<html lang="en">
	<head>
		<style type="text/css">
			body {font-family: arial;}
		</style>
		<meta http-equiv="content-type" content="text/html; charset=utf-8">
		<title>[Needs Title]</title>
	</head>
	<body>
		<form action="user.login_server.php" method="post">
			<label for="user">Username: </label><input type="text"     id="user" name="user"><br>
			<label for="pw"  >Password: </label><input type="password" id="pw"   name="pw"  ><br>
			<button type="submit">Login</button>
		</form>
		<br>
		If you don't have a user account, you may create one by clicking the "Register" button below.<br>
		<button type="button" onclick="window.location.href='register.php'">Register</button><br><br>
		<button type="button" onclick="window.location.href='index.php'">Back to Home</button>
	</body>
</html>
