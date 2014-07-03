function deleteHapmapConfirmation(hapmap,key){
	$("#h2_"+hapmap+"_delete").html("<b><font color=\"red\">[Are you sure?]</font><button type='button' onclick='deleteHapmap_yes(\""+hapmap+"\")'>Yes, delete.</button><button type='button' onclick='deleteHapmap_no(\""+hapmap+"\",\""+key+"\")'>No, cancel</button></b>");
	document.getElementById(key).style.display = 'none';
}

// Requires javascript value user to be defined globally
function deleteHapmap_yes(hapmap){
	$.ajax({
		url : 'php/hapmap.delete_server.php',
		type : 'post',
		data : {
			user: user,
			hapmap: hapmap
		},
		success : function(answer){
			if(answer == "COMPLETE"){
				window.location.href=window.location.href;
			}
		}
	});
}
function deleteHapmap_no(hapmap,key){
	$("#h2_"+hapmap+"_delete").html("");
	document.getElementById(key).style.display = 'inline';
}
