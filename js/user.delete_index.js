function deleteUserConfirmation(user,key){
	$("#u_"+user+"_delete").html("<b>This action will <font color=\"red\">PERMANENTLY</font> delete this user account and all associated files.<br>There are no backups to restore from.</br>Are you absolutely sure? <button type='button' onclick='deleteUser_yes(\""+user+"\")'>Yes, delete this user.</button><button type='button' onclick='deleteUser_no(\""+user+"\",\""+key+"\")'>No, don't delete this user.</button></b>");
	document.getElementById(key).style.display = 'none';
}

// Requires javascript value user to be defined globally
function deleteUser_yes(user){
	$.ajax({
		url : 'php/user.delete_server.php',
		type : 'post',
		data : {
			user : user
		},
		success : function(answer){
			if(answer == "COMPLETE"){
				window.location.href=window.location.href;
			}
		}
	});
}
function deleteUser_no(user,key){
	$("#u_"+user+"_delete").html("");
	document.getElementById(key).style.display = 'inline';
}
