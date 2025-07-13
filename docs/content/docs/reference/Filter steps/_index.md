---
bookCollapseSection: true
weight: 15
---

# Filter section

We can filter reads using the following filters.

Note the double `[[' in ``[filter]]``, since you're defining a list of filters.

Filters are interpreted eagerly - they're applied as early as possible to save on computation.

All filters have an `action` property that controls whether reads matching the filter
are removed ('remove') or the only ones kept ('keep').


{{<mynav>}}
