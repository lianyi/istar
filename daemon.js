#!/usr/bin/env node
var mongodb = require('mongodb');
new mongodb.MongoClient(new mongodb.Server(process.argv[2], 27017)).open(function(err, mongoClient) {
	if (err) throw err;
	var db = mongoClient.db('istar');
	db.authenticate(process.argv[3], process.argv[4], function(err, authenticated) {
		if (err) throw err;
		var idock = db.collection(process.argv[5]);
		if (process.argv[6] === undefined) {
			idock.count(function(err, count) {
				if (err) throw err;
				console.log(count);
				mongoClient.close(function(err, res) {
					if (err) throw err;
				});
			});
		} else {
			idock.find({}, {
				sort: {'submitted': 1},
				skip: process.argv[6],
			}).toArray(function(err, docs) {
				if (err) throw err;
				docs.forEach(function(doc) {
					console.log("%j", doc);
				});
				mongoClient.close(function(err, res) {
					if (err) throw err;
				});
			});
		}
	});
});
