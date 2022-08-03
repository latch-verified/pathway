import React from 'react';
import ReactDOM from 'react-dom/client';

import Select from 'react-select'
import createPanZoom from 'panzoom';
import EllipsisToolTip from "ellipsis-tooltip-react-chan";

import Plotly from 'plotly.js-dist';
import './index.css';

// Workflow code will inject JSON-stringified data into the report, using Jinja2
const injection = (name) => JSON.parse(name);
const PATHWAY_DATA = injection(`{{ pathwayData }}`)
const CONTRAST_DATA = injection(`{{ contrastData }}`)
const GENE_SETS_BY_PATHWAY_ID = injection(`{{ pathwayIdToGeneSets }}`)
const CORE_GENES = Array.from(new Set(PATHWAY_DATA.flatMap(x => x.coreEnrichedGenes)));
const PATHWAY_ID_TO_GENE_GROUPS = injection(`{{ pathwayIdToGeneGroups }}`)

const sorter = (a, b) => {
	const [af, bf] = [parseFloat(a), parseFloat(b)];
	if (isNaN(af) || isNaN(bf)) return String(a).localeCompare(String(b));
	else return af - bf;
}

const SortableTable = ({ headers, rows, rowIdKey, noneText, selected, selectable, onRowClick }) => {
	const [sort, setSort] = React.useState(null);

	function changeSort(key) {
		if (sort == null || sort[0] !== key) setSort([key, 'descending']);
		else setSort([key, sort[1] === 'ascending' ? 'descending' : 'ascending']);
	}

	let sortedRows = rows;
	if (sort != null) {
		if (sort[1] === 'ascending') {
			sortedRows.sort((a, b) => sorter(a[sort[0]], b[sort[0]]));
		} else {
			sortedRows.sort((a, b) => -sorter(a[sort[0]], b[sort[0]]));
		}
	}

	return (
		<div className="table-container">
			<table>
				<thead>
					<tr>
						{headers.map(({name, key}) => (
							<th key={key} onClick={() => changeSort(key)}>
								{name}
								{sort && sort[0] === key && sort[1] === 'ascending' && '(A)'}
								{sort && sort[0] === key && sort[1] === 'descending' && '(D)'}
							</th>
						))}
					</tr>
				</thead>
				<tbody>
					{sortedRows && sortedRows.length > 0 && sortedRows.map(row => {
						const a = selectable ? 'selectable' : '';
						const s = selected === row[rowIdKey] ? 'selected' : '';
						return (
							<tr 
								key={row[rowIdKey]}
								className={`${a} ${s}`}
								onClick={() => onRowClick(row, row[rowIdKey])}
							>
								{headers.map(({key, transformer}) => (
									<td key={key}>{transformer ? transformer(row[key]) : row[key]}</td>
								))}
							</tr>
						);
					})}
					{(!sortedRows || sortedRows.length === 0) && noneText}
				</tbody>
			</table>
		</div>
	);
};

const PathwayTable = ({ selectedPathwayId, setSelectedPathwayId }) => {
	const [filterGenes, setFilterGenes] = React.useState([]);

	const pathways = React.useMemo(() => {
		if (filterGenes == null)
			return;
		else if (filterGenes.length === 0)
			return PATHWAY_DATA;
		else
			return PATHWAY_DATA.filter(x => filterGenes.every(y => x.coreEnrichedGenes.includes(y)));
	}, [filterGenes]);

	return (
		<div className="f c s">
			<h2>Pathway Enrichment Table (KEGG only)</h2>
			<Select
				options={CORE_GENES.map(x => ({label: x, value: x}))}
				placeholder="Filter row by gene..."
				isMulti={true}
				onChange={values => setFilterGenes(values.map(x => x.value))}
			/>
			<br />
			<SortableTable
				headers={[
					{name: 'Pathway Name', key: 'pathwayName'},
					{name: 'p-Value', key: 'pValue'},
					{name: 'p-Adjusted', key: 'pAdjusted'},
					{name: 'ES', key: 'enrichmentScore'},
					{name: 'NES', key: 'normalizedEnrichmentScore'},
					{name: 'Gene Set Size', key: 'geneSetSize'},
					{name: 'Leading Edge Size', key: 'leadingEdgeSize'},
					{name: 'Core Enrichment', key: 'coreEnrichedGenes', transformer: x => x.join(" ")},
				]}
				rowIdKey="pathwayId"
				rows={pathways}
				noneText="No pathways found."
				selectable={true}
				selected={selectedPathwayId}
				onRowClick={(_, pathwayId) => setSelectedPathwayId(pathwayId)}
			/>
		</div>
	);
}

const pluralize = (value, n) => n === 1 ? value : `${value}s`;

const ViewGeneSelector = ({gene, selected, addViewGene, removeViewGene}) => (
	<input
		type="checkbox"
		checked={selected}
		onChange={(e) => {
			if (e.currentTarget.checked) {
				addViewGene(gene);
			} else {
				removeViewGene(gene);
			}
		}}
	/>
);

const GeneAnnotation = ({ core, view, genes, selected, addViewGene, removeViewGene }) => {
	const { x, y, width, height } = view;
	return (
		<>
			<div
				className={`map-area ${core ? 'map-area-core' : ''}`}
				style={{
					width: `${width}px`,
					height: `${height}px`,
					left: `${x - 2}px`,
					top: `${y - 2}px`,
				}}
			/>
			<div
				className="map-area-tag"
				style={{
					top: `${y}px`,
					left: `${x + width}px`,
				}}
			>
				<div className="select-names">
					Select {pluralize('gene', genes.length)}
				</div>
				{genes.map(({ name: gene }) => 
					<div>
						<ViewGeneSelector
							selected={selected.includes(gene)}
							gene={gene}
							addViewGene={addViewGene}
							removeViewGene={removeViewGene}
						/>
						{gene}
					</div>
				)}
			</div>
		</>
	);
}

const AnnotatedPathview = ({ pathwayId, selected, addViewGene, removeViewGene }) => {
	const ref = React.useRef();

	React.useEffect(() => {
		if (ref?.current == null) return;
		createPanZoom(ref.current);
	}, [ref]);

	return (
		<div class="f c s">
			<h2>Pathway Viewer (KEGG only)</h2>
			{/* <div id="pathview-container"> */}
			<h3 id="pathview-pathway-name">
				Select a pathway in the <b>Pathway Enrichment Table</b> to visualize a pathway graph.
			</h3>
			<h3 id="pathview-pathway-name">
				Hover over any gene to add or remove them from the <b>LFC plot</b> below.
			</h3>
			<div id="pathview-body">
				<div ref={ref} style={{position: 'relative'}}>
					{pathwayId && <img src={document.getElementById(`pathview-${pathwayId}`).getAttribute('src')} />}
					{pathwayId && PATHWAY_ID_TO_GENE_GROUPS[pathwayId].map(x =>
						<GeneAnnotation
							{...x}
							selected={selected}
							addViewGene={addViewGene}
							removeViewGene={removeViewGene}
						/>
					)}
				</div>
			</div>
		</div>
	);
}

const fetchGeneInfoByEntrezIds = async (ids) => {
	const body = {
		ids: ids.join(","),
		fields: ['geneid', 'name', 'type_of_gene', 'summary'].join(","),
	};
	const req = await fetch(`https://mygene.info/v3/gene`, {
		method: 'POST',
		headers: {'Content-Type': 'application/json'},
		body: JSON.stringify(body),
	});
	return await req.json();
}

const PathwayGeneTable = ({ selectedPathwayId, selectedGenes, addViewGene, removeViewGene }) => {
	const [filterGenes, setFilterGenes] = React.useState(null);
	const [store, setStore] = React.useState({});

	const genes = GENE_SETS_BY_PATHWAY_ID[selectedPathwayId];
	
	React.useEffect(() => {
		if (genes == null)
			return;

		const newIds = genes[0].filter(x => !store.hasOwnProperty(x));

		if (newIds.length === 0) 
			return;
		
		fetchGeneInfoByEntrezIds(newIds).then(data => {
			const newStore = store == null ? {} : {...store};
			for (const {name, query: id, summary, type_of_gene} of data) 
				newStore[id] = {name, type: type_of_gene, summary};
			setStore(newStore);
		});
	// eslint-disable-next-line react-hooks/exhaustive-deps
	}, [genes]);

	const rows = React.useMemo(() => {
		return genes?.[0].map((entrezId, i) => {
			const gene = genes[1][i];
			const contrastRow = CONTRAST_DATA[gene];
			const storeData = store[entrezId];
			const sd = (key) => {
				if (!storeData) return 'Loading';
				else if (storeData[key] == null) return 'Unavailable';
				else return storeData[key];
			}
			return {
				'name': gene,
				'lfc': contrastRow?.[0] ?? 'NA',
				'pValue': contrastRow?.[1] ?? 'NA',
				'pAdjusted': contrastRow?.[2] ?? 'NA',
				'type': sd('type'),
				'summary': sd('summary'),
			};
		}).filter(item => {
			return filterGenes == null || filterGenes.length === 0 || filterGenes.includes(item.name);
		});
	}, [filterGenes, genes, store]);

	const set = React.useMemo(() => new Set(selectedGenes), [selectedGenes]);

	return (
		<>
			<Select
				options={genes?.[1].map(x => ({label: x, value: x}))}
				placeholder="Search for gene..."
				isMulti={true}
				onChange={values => setFilterGenes(values.map(x => x.value))}
			/>
			<br />
			<SortableTable
				headers={[
					{name: 'Gene Name', key: 'name'},
					{name: 'LFC', key: 'lfc'},
					{name: 'p-val', key: 'pValue'},
					{name: 'p-adj', key: 'pAdjusted'},
					{name: 'Broad Class', key: 'type'},
					{name: 'Description', key: 'summary', transformer: item => (
						<div className="overflow-text-container">
							<EllipsisToolTip options={{place: 'top'}}>
								{item}
							</EllipsisToolTip>
						</div>
					)},
					{name: 'In LFC plot', key: 'name', transformer: item => (
						<ViewGeneSelector
							selected={set.has(item)}
							gene={item}
							addViewGene={addViewGene}
							removeViewGene={removeViewGene}
						/>
					)},
				]}
				rowIdKey="name"
				rows={rows}
				noneText="Select a pathway from the Pathway Enrichment Table to view its genes"
			/>
		</>
	);
}

const Plot = ({ data, layout }) => {
	const ref = React.useRef();

	React.useEffect(() => {
		Plotly.newPlot(ref.current, data, layout);
	}, [data, layout]);

	return <div ref={ref} />;
}

const LogFoldChangePlot = ({ selectedGenes }) => {
	const sortedGenes = React.useMemo(() => {
		return [...selectedGenes].sort((a, b) => -sorter(CONTRAST_DATA[b][0], CONTRAST_DATA[a][0]));
	}, [selectedGenes])
	
	const lfcs = sortedGenes.map(x => CONTRAST_DATA[x][0])

	return (
		<div className="plot-container">
			<Plot
				data={[
					{
						x: sortedGenes,
						y: lfcs,
						marker: {
							color: lfcs.map(x => x > 0 ? 'red' : 'blue'),
						},
						type: 'bar',
					}
				]}
				layout={{
					title: 'Log fold change (LFC)'
				}}
			/>
		</div>
	);
}

const Report = () => {
	const [selectedPathwayId, setSelectedPathwayId] = React.useState(null);
	const [viewGenes, setViewGenes] = React.useState([]);

	const addViewGene = (gene) => {
		setViewGenes([...viewGenes, gene]);
	}

	const removeViewGene = (gene) => {
		const index = viewGenes.indexOf(gene);
		setViewGenes([...viewGenes.slice(0, index), ...viewGenes.slice(index + 1)]);
	}

	React.useEffect(() => {
		if (selectedPathwayId == null) return;
		const coreGenes = PATHWAY_ID_TO_GENE_GROUPS[selectedPathwayId].flatMap(
			group => group.genes.filter(x => x.core).map(x => x.name)
		);
		setViewGenes([...new Set(coreGenes)]);
	}, [selectedPathwayId])

	return (
		<>
			<PathwayTable
			 	selectedPathwayId={selectedPathwayId}
				setSelectedPathwayId={setSelectedPathwayId}
			/>
			<AnnotatedPathview
				pathwayId={selectedPathwayId}
				selected={viewGenes}
				addViewGene={addViewGene}
				removeViewGene={removeViewGene}
			/>
			<div className="f c s">
				<h2>Log fold change (LFC) plot</h2>
				<h3 id="lfc-help">
					Select genes in the pathway viewer or the pathway gene set table to render the LFC plot.
				</h3>
				<LogFoldChangePlot selectedGenes={viewGenes} />
				<br />
				<button onClick={() => setViewGenes([])}>Clear all genes</button>
			</div>
			<div className="f c s">
				<h2>Pathway Gene Set Table</h2>
				<PathwayGeneTable
					selectedPathwayId={selectedPathwayId}
					selectedGenes={viewGenes}
					addViewGene={addViewGene}
					removeViewGene={removeViewGene}
				/>
			</div>
		</>
	);
}

const root = ReactDOM.createRoot(document.getElementById('react-entrypoint'));
root.render(<React.StrictMode><Report /></React.StrictMode>);
