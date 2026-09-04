document.addEventListener("DOMContentLoaded", function() {

  // --- CONFIG ---
  const REMOVE_UNNUMBERED_FROM_TOC = true; 

  // --- HELPER: FORMAT NUMBERS ---
  function formatLinkText(link) {
    let container = link.querySelector('.menu-text') || link;
    if (container.querySelector('.toc-num')) return false;

    const text = container.innerText;
    if (!text) return false;

    const match = text.match(/^([\d\w\.]+\.?)(\s+)(.*)/);

    if (match) {
      container.innerHTML = `<span class="toc-num">${match[1]}</span><span class="toc-text">${match[3]}</span>`;
      link.classList.add('is-numbered');
      return true;
    } else {
      if (REMOVE_UNNUMBERED_FROM_TOC && link.closest('.injected-page-toc')) {
          const li = link.closest('li');
          if (li) li.remove();
          return false; 
      }
      link.classList.add('is-unnumbered');
      return true;
    }
  }

  // --- PART A: FORMAT CHAPTERS ---
  function formatChapters() {
    const items = document.querySelectorAll('#quarto-sidebar .sidebar-item');
    items.forEach(li => {
      const link = li.querySelector('.sidebar-link');
      const toggleBtn = li.querySelector('.sidebar-item-toggle');

      if (link) {
        formatLinkText(link);
        
        if (toggleBtn) {
          link.classList.add('has-children');
          if (li.classList.contains('sidebar-item-collapsed')) {
             link.classList.add('collapsed');
          } else {
             link.classList.add('expanded');
          }

          if (!link.dataset.hasListener) {
            link.addEventListener('click', (e) => {
               e.stopPropagation(); 
               if (link.classList.contains('collapsed')) {
                   link.classList.remove('collapsed');
                   link.classList.add('expanded');
               } else {
                   link.classList.add('collapsed');
                   link.classList.remove('expanded');
               }
               toggleBtn.click();
            });
            link.dataset.hasListener = "true";
          }
        }
      }
    });
  }

  // --- PART B: INJECT PAGE TOC ---
  function moveToc() {
    const toc = document.querySelector('nav[role="doc-toc"]');
    const activeLink = document.querySelector('.sidebar-link.active');

    if (toc && activeLink) {
      const parentLi = activeLink.closest('.sidebar-item');
      if (parentLi && !parentLi.querySelector('.injected-page-toc')) {
        const tocList = toc.querySelector('ul');
        if (tocList) {
            const clonedList = tocList.cloneNode(true);
            clonedList.className = 'injected-page-toc';

            const links = clonedList.querySelectorAll('a');
            links.forEach(link => formatLinkText(link));

            const parents = clonedList.querySelectorAll('li');
            parents.forEach(li => {
                if (li.querySelector('ul')) { 
                    const link = li.querySelector('a');
                    const subList = li.querySelector('ul');
                    if (link && subList) {
                        link.classList.add('has-children');
                        link.classList.add('collapsed'); 
                        subList.style.display = 'none';

                        if (!link.dataset.hasListener) {
                            link.addEventListener('click', (e) => {
                                e.stopPropagation();
                                if (subList.style.display === 'none') {
                                    subList.style.display = 'block';
                                    link.classList.remove('collapsed');
                                    link.classList.add('expanded');
                                } else {
                                    subList.style.display = 'none';
                                    link.classList.remove('expanded');
                                    link.classList.add('collapsed');
                                }
                            });
                            link.dataset.hasListener = "true";
                        }
                    }
                }
            });
            parentLi.appendChild(clonedList);
            syncActiveState(clonedList); 
        }
      }
    }
  }

  // --- PART C: LIVE SYNC (Title Priority Fix) ---
  function syncActiveState(injectedToc) {
    const tocLinks = injectedToc.querySelectorAll('a');
    tocLinks.forEach(l => l.classList.remove('active')); 
    tocLinks.forEach(l => l.classList.remove('toc-active')); 

    const headers = document.querySelectorAll('main h1, main h2, main h3, main h4, main h5, main h6');
    const idToLink = {};
    
    tocLinks.forEach(link => {
        const href = link.getAttribute('href');
        if(href && href.startsWith('#')) {
            const id = href.substring(1);
            idToLink[id] = link;
            idToLink[decodeURIComponent(id)] = link; 
        }
    });

    const titleObserver = new IntersectionObserver((entries) => {
        entries.forEach(entry => {
            if (entry.isIntersecting) {
                tocLinks.forEach(l => l.classList.remove('toc-active'));
            }
        });
    }, { rootMargin: '0px 0px -50% 0px' }); 

    const mainTitle = document.querySelector('main h1');
    if (mainTitle) titleObserver.observe(mainTitle);

    const contentObserverOptions = {
        root: null,
        rootMargin: '0px 0px -80% 0px', 
        threshold: 0
    };

    const contentObserver = new IntersectionObserver((entries) => {
        entries.forEach(entry => {
            if (entry.isIntersecting) {
                if (window.scrollY < 100) return; 

                const id = entry.target.getAttribute('id');
                const activeSubLink = idToLink[id];
                
                if (activeSubLink) {
                    tocLinks.forEach(l => l.classList.remove('toc-active'));
                    activeSubLink.classList.add('toc-active');
                    
                    const parentUl = activeSubLink.closest('ul');
                    if (parentUl && parentUl.style.display === 'none') {
                           parentUl.style.display = 'block';
                           const pLink = parentUl.parentElement.querySelector('a.has-children');
                           if (pLink) {
                               pLink.classList.remove('collapsed');
                               pLink.classList.add('expanded');
                           }
                    }
                }
            }
        });
    }, contentObserverOptions);

    headers.forEach(header => {
        if (header.tagName !== 'H1') {
            contentObserver.observe(header);
        }
    });
  }

  // --- PART D: SIDEBAR TOGGLE ---
  function createSidebarToggle() {
    if (document.getElementById('custom-sidebar-toggle')) return;
    const btn = document.createElement("button");
    btn.id = "custom-sidebar-toggle";
    btn.innerHTML = '&#9776;'; 
    btn.title = "Toggle Sidebar";
    document.body.insertAdjacentElement('afterbegin', btn);
    btn.addEventListener("click", function() {
      document.body.classList.toggle("sidebar-closed");
    });
  }

  // --- PART E: GITHUB SOURCE LINK ---
  function addGithubLink() {
    if (document.getElementById('injected-github-link')) return;

    // Look for Quarto's natively generated GitHub links (either in toc-actions or navbar tools)
    const githubNode = document.querySelector('a.toc-action[href*="github.com"], a.quarto-navigation-tool[href*="github.com"]');
    
    if (githubNode) {
      const githubUrl = githubNode.getAttribute('href');
      const activeLink = document.querySelector('.sidebar-link.active');
      const parentLi = activeLink ? activeLink.closest('.sidebar-item') : null;

      if (parentLi) {
        const linkContainer = document.createElement('div');
        linkContainer.id = 'injected-github-link';
        // Inline styles to match typical sidebar formatting
        linkContainer.style.cssText = 'margin-top: 1rem; padding-top: 0.8rem; padding-left: 1rem; border-top: 1px solid var(--bs-border-color, rgba(0,0,0,0.1)); font-size: 0.85em; opacity: 0.85;';

        linkContainer.innerHTML = `
          <a href="${githubUrl}" target="_blank" rel="noopener noreferrer" style="text-decoration: none; color: inherit; display: flex; align-items: center; gap: 6px; transition: opacity 0.2s;" onmouseover="this.style.opacity='1'" onmouseout="this.style.opacity='0.85'">
            <svg style="width: 16px; height: 16px;" fill="currentColor" viewBox="0 0 24 24" aria-hidden="true">
              <path fill-rule="evenodd" d="M12 2C6.477 2 2 6.484 2 12.017c0 4.425 2.865 8.18 6.839 9.504.5.092.682-.217.682-.483 0-.237-.008-.868-.013-1.703-2.782.605-3.369-1.343-3.369-1.343-.454-1.158-1.11-1.466-1.11-1.466-.908-.62.069-.608.069-.608 1.003.07 1.531 1.032 1.531 1.032.892 1.53 2.341 1.088 2.91.832.092-.647.35-1.088.636-1.338-2.22-.253-4.555-1.113-4.555-4.951 0-1.093.39-1.988 1.029-2.688-.103-.253-.446-1.272.098-2.65 0 0 .84-.27 2.75 1.026A9.564 9.564 0 0112 6.844c.85.004 1.705.115 2.504.337 1.909-1.296 2.747-1.027 2.747-1.027.546 1.379.202 2.398.1 2.651.64.7 1.028 1.595 1.028 2.688 0 3.848-2.339 4.695-4.566 4.943.359.309.678.92.678 1.855 0 1.338-.012 2.419-.012 2.747 0 .268.18.58.688.482A10.019 10.019 0 0022 12.017C22 6.484 17.522 2 12 2z" clip-rule="evenodd"></path>
            </svg>
            <span>View source on GitHub</span>
          </a>
        `;
        
        parentLi.appendChild(linkContainer);
      }
    }
  }

  // --- EXECUTION ---
  function run() {
    try { createSidebarToggle(); } catch(e) { console.error(e); }
    try { formatChapters(); } catch(e) { console.error(e); }
    try { moveToc(); } catch(e) { console.error(e); }
    try { addGithubLink(); } catch(e) { console.error(e); }
    document.querySelector('#quarto-sidebar .sidebar-menu-container')?.classList.add('loaded');
  }

  run();
  setTimeout(run, 500); 
});
