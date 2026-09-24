/* pkgdown loads this file in <head>, so wait for the page before touching it */
document.addEventListener('DOMContentLoaded', function () {

  /* Back-to-top button */
  var top = document.createElement('button');
  top.id = 'back-to-top';
  top.title = 'Back to top';
  top.innerHTML = '&#8679;'; /* up arrow */
  document.body.appendChild(top);
  window.addEventListener('scroll', function () {
    top.style.display = window.scrollY > 300 ? 'block' : 'none';
  });
  top.addEventListener('click', function () {
    window.scrollTo({ top: 0, behavior: 'smooth' });
  });

  /* Lab bar colour toggle: stores the same "lab-theme" key as avisrilab.org,
     so the choice carries between the lab site and these docs */
  var btn = document.getElementById('labTheme');
  if (!btn) return;
  var html = document.documentElement;
  var sun = '<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="5"/><line x1="12" y1="1" x2="12" y2="3"/><line x1="12" y1="21" x2="12" y2="23"/><line x1="4.22" y1="4.22" x2="5.64" y2="5.64"/><line x1="18.36" y1="18.36" x2="19.78" y2="19.78"/><line x1="1" y1="12" x2="3" y2="12"/><line x1="21" y1="12" x2="23" y2="12"/><line x1="4.22" y1="19.78" x2="5.64" y2="18.36"/><line x1="18.36" y1="5.64" x2="19.78" y2="4.22"/></svg>';
  var moon = '<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"><path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"/></svg>';
  function icon() { btn.innerHTML = html.getAttribute('data-bs-theme') === 'light' ? moon : sun; }
  icon();
  btn.addEventListener('click', function () {
    var next = html.getAttribute('data-bs-theme') === 'light' ? 'dark' : 'light';
    html.setAttribute('data-bs-theme', next);
    try { localStorage.setItem('lab-theme', next); } catch (e) {}
    icon();
  });
});
