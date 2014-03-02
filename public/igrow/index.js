$(function() {

	// Initialize pager
	var pager = $('#pager');
	pager.pager('init', [ 'Description', 'Submitted on', 'Status', 'Progress', 'Result' ], function(job) {
		var status, progress, result = '<a href="iview/?' + job._id + '"><img src="/iview/logo.png" alt="iview"></a>';
		if (!job.scheduled) {
			status = 'Queued for execution';
			progress = 0;
		} else if (!job.done) {
			status = 'Execution in progress';
			progress = job.generation / job.num_generations;
		} else {
			status = 'Done on ' + $.format.date(new Date(job.done), 'yyyy/MM/dd HH:mm:ss');
			progress = 1;
			result += '<a href="jobs/' + job._id + '/log.csv.gz"><img src="/excel.png" alt="log.csv.gz"></a><a href="jobs/' + job._id + '/ligands.pdbqt.gz"><img src="/molecule.png" alt="ligands.pdbqt.gz"></a>';
		}
		return [
			job.description,
			$.format.date(new Date(job.submitted), 'yyyy/MM/dd HH:mm:ss'),
			status,
			(100 * progress).toFixed(5) + '%',
			result
		];
	});

	// Refresh the table of jobs and its pager every second
	var jobs = [], skip = 0;
	var tick = function() {
		$.get('jobs', { skip: skip, count: jobs.length }, function(res) {
			if (res.length) {
				for (var i = skip; i < jobs.length; ++i) {
					var job = res[i - skip];
					jobs[i].done = job.done;
				}
				pager.pager('refresh', skip, jobs.length, 2, 5, false);
				if (res.length > jobs.length - skip) {
					var len = jobs.length;
					jobs = jobs.concat(res.slice(jobs.length - skip));
					pager.pager('source', jobs);
					pager.pager('refresh', len, jobs.length, 0, 5, true);
				}
				for (skip = jobs.length; skip && !jobs[skip - 1].done; --skip);
			}
			setTimeout(tick, 1000);
		});
	};
	tick();

	// Initialize tooltips
	$('.form-group a').tooltip();

	// Process submission
	var submit = $('#submit');
	submit.click(function() {
		// Hide tooltips
		$('.form-group a').tooltip('hide');
		// Do client side validation
//		var idock_id = $('#description').val();
		var generations = parseInt($('#generations').val());
		var email = $('#email').val();
		var v = new validator({
			email: email,
			description: description,
			generations: generations,
		});
		if (v
			.field('generations').message('must be an integer within [1, 8]').int().min(1).max(8)
			.field('email').message('must be valid').email()
			.failed()) {
			var keys = Object.keys(v.err);
			keys.forEach(function(key) {
				$('#' + key + '_label').tooltip('show');
			});
			$('#' + keys[0]).focus();
			return;
		}
		// Disable the submit button for a while
		submit.prop('disabled', true);
		// Post a new job with server side validation
		$.post('jobs', {
			generations: generations,
			description: description,
			email: email,
		}, function(res) {
			var keys = Object.keys(res);
			// If server side validation fails, show the tooltips
			if (keys.length) {
				keys.forEach(function(key) {
					$('#' + key + '_label').tooltip('show');
				});
			} else {
				$('html, body').animate({ scrollTop: pager.offset().top });
//				window.scrollTo(pager.offset().left, pager.offset().top);
			}
		}, 'json').always(function() {
			submit.prop('disabled', false);
		});
	});

	// Apply accordion to tutorials
	$('.ui-accordion').accordion({
		collapsible: true,
		active: false,
		heightStyle: 'content',
		activate: function(event, ui) {
			$('img', this).trigger('expand');
		}
	});
	$('.ui-accordion img').lazyload({
		event: 'expand',
		effect: 'fadeIn',
	});
});
