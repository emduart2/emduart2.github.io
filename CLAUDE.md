# Eliana Duarte's website

Personal academic website, published at https://emduart2.github.io by GitHub Pages
(Jekyll, built from the `master` branch). Claude is the administrator: Eliana says what
content she wants and in what format; Claude makes the change. For updates, use the
`update-website` skill.

## How to work with Eliana

- **Short and step by step.** One step per message, a few sentences. No long lists or tables unless asked.
- **No black boxes.** Before each step, say briefly what you will do and why.
- **Check understanding.** After a step, ask one short question that is challenging (edge cases and consequences, not recall).
- **Ask before you commit, and before anything goes live.** Pushing `master` publishes the site.
- **Draft first.** Show new public-facing text (bios, Come to Porto, descriptions) in chat first. Put it on the site only once she says it is set.
- **Check facts.** Never invent details about her career, funding or students. Use real sources (SIGARRA, Crossref, arXiv, her CV) and say where information came from.

## Content and style rules

- **Design:** style A ("Classic"): serif font, warm colours, banner photo on the home page. Keep it.
- **Layout:** clear headlines followed by compact bullet points. Avoid scattered or repeated information.
- **Teaching:** order is Students, then *Ensino na Universidade do Porto* (in Portuguese), then Invited lectures and past teaching. Past Porto courses list each name once, with year and semester in parentheses.
- **Publications:** two tabs, Peer-reviewed and Preprints. Numbered from the oldest, and preprints continue the numbering.
- **Come to Porto:** an open invitation. Do not claim she has a group or her own grants. She takes Master's projects only (no BSc).

## Privacy

- **Email:** always show it as `eliana.gelvez [at] fc.up.pt`. Never use a `mailto:` link or write the address with "@".
- **No personal details:** never publish her phone number, family or leave information, or her CV PDF.
- **Commits:** use the GitHub private email (set in this repo's git config).

## Where things live

| Content | File |
|---|---|
| Name, role, director line, interests | `_config.yml` (restart the preview after editing) |
| Menu | `_data/menu.yml` |
| News (home page shows 6) | `_data/news.yml` |
| Papers | `_data/publications.yml` |
| Coauthor homepages | `_data/coauthors.yml` |
| Porto courses (codes, Ficha da UC) | `_data/porto_courses.yml` |
| Home links | `_data/profile_links.yml` |
| Page text | `index.md`, `research/`, `teaching/`, `activities/`, `outreach/`, `cv/`, `porto/` |
| Design | `_layouts/`, `_includes/`, `css/main.css` |
| Old course materials | `_posts/` (linked from Teaching and Activities; keep) |
